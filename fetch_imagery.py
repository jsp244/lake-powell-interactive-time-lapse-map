"""Fetch yearly Landsat true-color frames of Lake Powell from Microsoft Planetary Computer.

For each year since 1984, finds the clearest summer scene(s) covering the lake
(mosaicking adjacent WRS-2 path/row scenes when needed), renders a true-color
JPEG cropped to the lake, and estimates water surface area via NDWI.

Writes frames/<year>.jpg, manifest.json, and manifest.js (for file:// viewing).
Resumable: skips years whose frame already exists unless --force.
"""

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import planetary_computer
import rasterio
from PIL import Image
from pystac_client import Client
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform_bounds
from shapely.geometry import box, shape
from shapely.ops import unary_union

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "landsat-c2-l2"
# Pre-1984: Landsat 1-3 MSS, Level-1 only (no blue band, no SWIR, 60 m)
MSS_COLLECTION = "landsat-c2-l1"
MSS_LAST_YEAR = 1983
MSS_BANDS = ["red", "green", "nir09"]
# Mean exoatmospheric irradiance (W/m2/um) per MSS band, for TOA reflectance
MSS_ESUN = {"green": 1848.0, "red": 1588.0, "nir09": 880.0}
# MSS detector-dropout streaks fool NDWI; full pool is ~653 km2, so anything
# above this is a broken scene -> try the next candidate date
MSS_MAX_PLAUSIBLE_KM2 = 700.0

# Lake Powell AOI (lon/lat): Glen Canyon Dam up to Hite, incl. San Juan/Escalante arms
AOI = (-111.62, 36.85, -110.35, 37.95)
DST_CRS = "EPSG:32612"  # UTM 12N
RES = 75.0  # meters/pixel output grid

# Fixed stretch so brightness is consistent across years (no flicker)
VMIN, VMAX, GAMMA = 0.012, 0.33, 1.55
# Water = MNDWI (green vs SWIR, haze-robust) positive, on QA-clear pixels
MNDWI_THRESHOLD = 0.0
# qa_pixel bits: 1 dilated cloud, 2 cirrus, 3 cloud, 4 cloud shadow
QA_BAD_BITS = (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4)

GDAL_ENV = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    GDAL_HTTP_MAX_RETRY="6",
    GDAL_HTTP_RETRY_DELAY="2",
    CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.TIF",
    GDAL_HTTP_MULTIPLEX="YES",
)

BANDS = ["red", "green", "blue", "swir16"]


def dst_grid():
    left, bottom, right, top = transform_bounds("EPSG:4326", DST_CRS, *AOI)
    width = int(round((right - left) / RES))
    height = int(round((top - bottom) / RES))
    transform = rasterio.transform.from_origin(left, top, RES, RES)
    return transform, width, height


def search_year(catalog, year, season, cloud_lt, platforms=None, collection=COLLECTION):
    query = {"eo:cloud_cover": {"lt": cloud_lt}}
    if platforms:
        query["platform"] = {"in": platforms}
    search = catalog.search(
        collections=[collection],
        bbox=AOI,
        datetime=f"{year}-{season[0]}/{year}-{season[1]}",
        query=query,
        limit=200,
    )
    return list(search.items())


def score_groups(items, aoi_poly):
    """Group items by (date, platform); return list of (cov, -max_cloud, key, group)
    sorted best-first."""
    groups = defaultdict(list)
    for it in items:
        d = it.datetime.date().isoformat()
        groups[(d, it.properties.get("platform", "?"))].append(it)
    scored = []
    for key, group in groups.items():
        union = unary_union([shape(it.geometry) for it in group])
        cov = union.intersection(aoi_poly).area / aoi_poly.area
        max_cloud = max(it.properties.get("eo:cloud_cover", 100) for it in group)
        scored.append((cov, -max_cloud, key, group))
    scored.sort(key=lambda t: (t[0], t[1]), reverse=True)
    return scored


def build_set(primary_idx, scored, aoi_poly):
    """Start from one group, fill AOI gaps with other groups nearest in time.
    Returns (chosen_items, coverage_ratio)."""
    best_cov, _, best_key, chosen = scored[primary_idx]
    chosen = list(chosen)
    covered = unary_union([shape(it.geometry) for it in chosen])
    cov = best_cov
    if cov < 0.995:
        best_date = best_key[0]
        rest = [s for i, s in enumerate(scored) if i != primary_idx]
        rest.sort(
            key=lambda t: abs(
                (np.datetime64(t[2][0]) - np.datetime64(best_date)).astype(int)
            )
        )
        for _, _, _, group in rest:
            union = unary_union([shape(it.geometry) for it in group])
            gain = union.difference(covered).intersection(aoi_poly).area / aoi_poly.area
            if gain > 0.005:
                chosen = chosen + group
                covered = unary_union([covered, union])
                cov = covered.intersection(aoi_poly).area / aoi_poly.area
            if cov >= 0.995 or len(chosen) >= 5:
                break
    return chosen, cov


def read_qa_bad_frac(items, transform, width, height):
    """Mosaic qa_pixel at 1/4 resolution over the AOI; return fraction of AOI
    pixels flagged cloud/cirrus/shadow."""
    w, h = width // 4, height // 4
    t = transform * transform.scale(4, 4)
    qa = np.zeros((h, w), dtype=np.uint16)
    for it in items:
        signed = planetary_computer.sign(it)
        with rasterio.open(signed.assets["qa_pixel"].href) as src:
            with WarpedVRT(
                src, crs=DST_CRS, transform=t, width=w, height=h,
                resampling=Resampling.nearest, src_nodata=1, nodata=1,
            ) as vrt:
                data = vrt.read(1)
        hole = qa <= 1
        qa[hole] = data[hole]
    valid = qa > 1
    if not valid.any():
        return 1.0
    bad = (qa & QA_BAD_BITS) > 0
    return float((bad & valid).sum() / valid.sum())


def pick_items(items, aoi_poly, transform, width, height):
    """Evaluate up to 4 candidate scene-sets; score each by cloud/shadow inside
    the AOI (per the QA band). Returns ranked [(aoi_bad_frac, -coverage, items)],
    best first."""
    if not items:
        return []
    scored = score_groups(items, aoi_poly)
    candidates, seen = [], set()
    for i in range(min(4, len(scored))):
        chosen, cov = build_set(i, scored, aoi_poly)
        sig = tuple(sorted(it.id for it in chosen))
        if sig not in seen:
            seen.add(sig)
            candidates.append((chosen, cov))
    full = [c for c in candidates if c[1] >= 0.985] or candidates
    results = []
    for chosen, cov in full:
        try:
            bad = read_qa_bad_frac(chosen, transform, width, height)
        except Exception:
            bad = max(
                it.properties.get("eo:cloud_cover", 100) for it in chosen
            ) / 100.0
        results.append((bad, -cov, chosen))
        if bad < 0.01:
            break  # essentially clear over the lake; stop probing
    results.sort(key=lambda t: (t[0], t[1]))
    return results


def read_mosaic(items, transform, width, height):
    """Return (bands[4,H,W] float32 surface reflectance, valid mask, clear mask).
    Items are mosaicked in order: first-valid wins. clear = valid and not flagged
    cloud/cirrus/shadow by the QA band."""
    dn = np.zeros((len(BANDS), height, width), dtype=np.uint16)
    qa = np.zeros((height, width), dtype=np.uint16)
    for it in items:
        signed = planetary_computer.sign(it)
        for bi, band in enumerate(BANDS):
            href = signed.assets[band].href
            with rasterio.open(href) as src:
                with WarpedVRT(
                    src,
                    crs=DST_CRS,
                    transform=transform,
                    width=width,
                    height=height,
                    resampling=Resampling.bilinear,
                    src_nodata=0,
                    nodata=0,
                ) as vrt:
                    data = vrt.read(1)
            hole = dn[bi] == 0
            dn[bi][hole] = data[hole]
        with rasterio.open(signed.assets["qa_pixel"].href) as src:
            with WarpedVRT(
                src, crs=DST_CRS, transform=transform, width=width, height=height,
                resampling=Resampling.nearest, src_nodata=1, nodata=1,
            ) as vrt:
                qdata = vrt.read(1)
        qhole = qa <= 1
        qa[qhole] = qdata[qhole]
    valid = np.all(dn > 0, axis=0)
    clear = valid & ((qa & QA_BAD_BITS) == 0)
    sr = dn.astype(np.float32) * 0.0000275 - 0.2
    sr = np.clip(sr, 0.0, 1.0)
    return sr, valid, clear


def read_mosaic_mss(items, transform, width, height):
    """MSS path: per scene, convert uint8 DN -> TOA reflectance using the scene's
    radiometric calibration and sun elevation, then mosaic (first-valid wins).
    Returns (bands[3,H,W]: red, green, nir09; valid; clear)."""
    refl = np.zeros((len(MSS_BANDS), height, width), dtype=np.float32)
    filled = np.zeros((height, width), dtype=bool)
    qa = np.zeros((height, width), dtype=np.uint16)
    for it in items:
        signed = planetary_computer.sign(it)
        sun_el = np.radians(it.properties.get("view:sun_elevation", 60.0))
        item_valid = None
        item_refl = np.zeros_like(refl)
        for bi, band in enumerate(MSS_BANDS):
            asset = signed.assets[band]
            rb = (asset.extra_fields.get("raster:bands") or [{}])[0]
            scale = rb.get("scale", 1.0)
            offset = rb.get("offset", 0.0)
            with rasterio.open(asset.href) as src:
                with WarpedVRT(
                    src, crs=DST_CRS, transform=transform, width=width,
                    height=height, resampling=Resampling.bilinear,
                    src_nodata=0, nodata=0,
                ) as vrt:
                    dn = vrt.read(1)
            radiance = dn.astype(np.float32) * scale + offset
            item_refl[bi] = (
                np.pi * radiance / (MSS_ESUN[band] * np.sin(sun_el))
            )
            bv = dn > 0
            item_valid = bv if item_valid is None else (item_valid & bv)
        with rasterio.open(signed.assets["qa_pixel"].href) as src:
            with WarpedVRT(
                src, crs=DST_CRS, transform=transform, width=width, height=height,
                resampling=Resampling.nearest, src_nodata=1, nodata=1,
            ) as vrt:
                qdata = vrt.read(1)
        hole = ~filled & item_valid
        refl[:, hole] = item_refl[:, hole]
        filled |= item_valid
        qhole = qa <= 1
        qa[qhole] = qdata[qhole]
    refl = np.clip(refl, 0.0, 1.0)
    clear = filled & ((qa & QA_BAD_BITS) == 0)
    return refl, filled, clear


def render(sr, valid, mss=False):
    if mss:
        # Render-only dark-object subtraction (anchored on the lake) so hazy
        # MSS TOA roughly matches the atmospherically-corrected L2 era, then
        # simulated natural color (no blue band: synthesize from green).
        red, green = sr[0].copy(), sr[1].copy()
        if valid.any():
            for band in (red, green):
                p1 = float(np.percentile(band[valid], 1))
                band -= max(0.0, p1 - 0.01)
        np.clip(red, 0.0, 1.0, out=red)
        np.clip(green, 0.0, 1.0, out=green)
        blue = green * 0.72
        rgb = np.stack([red, green, blue])
    else:
        rgb = sr[:3]
    x = np.clip((rgb - VMIN) / (VMAX - VMIN), 0.0, 1.0)
    x = x ** (1.0 / GAMMA)
    img = (x * 255).astype(np.uint8)
    img[:, ~valid] = 0
    return np.moveaxis(img, 0, -1)


def water_area_km2(sr, clear, mss=False):
    # TM/OLI: MNDWI (green vs SWIR). MSS has no SWIR: NDWI (green vs NIR).
    green = sr[1]
    other = sr[2] if mss else sr[3]
    denom = green + other
    with np.errstate(divide="ignore", invalid="ignore"):
        mndwi = np.where(denom > 0, (green - other) / denom, -1)
    water = (mndwi > MNDWI_THRESHOLD) & clear
    px_km2 = (RES * RES) / 1e6
    area = float(water.sum()) * px_km2
    clear_frac = float(clear.mean())
    corrected = area / clear_frac if 0.5 < clear_frac < 0.98 else area
    return round(corrected, 1), round(clear_frac, 3)


def process_year(catalog, year, aoi_poly, transform, width, height):
    mss = year <= MSS_LAST_YEAR
    if mss:
        attempts = [
            (("06-01", "10-15"), 10, None),
            (("06-01", "10-15"), 30, None),
            (("05-01", "10-31"), 60, None),
            (("01-01", "12-31"), 80, None),  # sparse years (1982): take anything
        ]
        collection = MSS_COLLECTION
    else:
        natural = ["landsat-5", "landsat-8", "landsat-9"]
        attempts = [
            (("06-01", "10-15"), 10, natural),
            (("06-01", "10-15"), 25, natural),
            (("05-01", "10-31"), 40, natural),
            (("06-01", "10-15"), 25, None),  # allow Landsat 7 (SLC-off stripes)
            (("05-01", "10-31"), 60, None),
        ]
        collection = COLLECTION
    all_items, seen_items, ranked = [], set(), []
    for season, cloud, platforms in attempts:
        for it in search_year(catalog, year, season, cloud, platforms, collection):
            if it.id not in seen_items:
                seen_items.add(it.id)
                all_items.append(it)
        ranked = [
            r
            for r in pick_items(all_items, aoi_poly, transform, width, height)
            if -r[1] >= 0.7
        ]
        if ranked and ranked[0][0] < 0.03 and -ranked[0][1] >= 0.985:
            break
    if not ranked:
        return None
    # MSS detector-dropout streaks can masquerade as water; if the area is
    # physically impossible for the era, retry with the tried dates excluded.
    # (Pool first neared full in 1979; before that ~600 km2 is the ceiling.)
    max_area = (600.0 if year <= 1978 else MSS_MAX_PLAUSIBLE_KM2) if mss else None
    tries, excluded = [], set()
    for _ in range(3 if mss else 1):
        pool = [it for it in all_items if it.id not in excluded]
        cands = [
            r
            for r in pick_items(pool, aoi_poly, transform, width, height)
            if -r[1] >= 0.7
        ]
        if not cands:
            break
        aoi_bad, neg_cov, chosen = cands[0]
        if mss:
            sr, valid, clear = read_mosaic_mss(chosen, transform, width, height)
        else:
            sr, valid, clear = read_mosaic(chosen, transform, width, height)
        img = render(sr, valid, mss=mss)
        area, clear_frac = water_area_km2(sr, clear, mss=mss)
        tries.append((area, img, clear_frac, aoi_bad, neg_cov, chosen))
        if max_area is None or area <= max_area:
            break
        excluded.update(it.id for it in chosen)
    # accept the first plausible try, else the least-inflated one
    if max_area is not None:
        tries.sort(key=lambda t: (t[0] > max_area, t[0]))
    area, img, clear_frac, aoi_bad, neg_cov, chosen = tries[0]
    cov = -neg_cov
    dates = sorted({it.datetime.date().isoformat() for it in chosen})
    platforms_used = sorted({it.properties.get("platform", "?") for it in chosen})
    clouds = [it.properties.get("eo:cloud_cover", -1) for it in chosen]
    return {
        "img": img,
        "meta": {
            "year": year,
            "dates": dates,
            "platforms": platforms_used,
            "sensor": "MSS" if mss else "TM/OLI",
            "cloud": round(max(clouds), 1),
            "aoi_cloud": round(aoi_bad, 3),
            "coverage": round(cov, 3),
            "clear_frac": clear_frac,
            "area_km2": area,
            "scenes": len(chosen),
        },
    }


def write_manifest(out_dir, entries):
    manifest = {
        "aoi": AOI,
        "res_m": RES,
        "years": sorted(entries, key=lambda e: e["year"]),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=1))
    (out_dir / "manifest.js").write_text(
        "const TIMELAPSE = " + json.dumps(manifest) + ";\n"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", default="1972-2026")
    ap.add_argument("--only", type=int)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.only:
        years = [args.only]
    else:
        a, b = args.years.split("-")
        years = list(range(int(a), int(b) + 1))

    out_dir = Path(__file__).parent
    frames = out_dir / "frames"
    frames.mkdir(exist_ok=True)
    aoi_poly = box(*AOI)
    transform, width, height = dst_grid()
    print(f"Output grid: {width}x{height} px at {RES} m", flush=True)

    catalog = Client.open(STAC_URL)
    entries = []
    manifest_path = out_dir / "manifest.json"
    if manifest_path.exists():
        entries = json.loads(manifest_path.read_text())["years"]

    with rasterio.Env(**GDAL_ENV):
        for year in years:
            frame_path = frames / f"{year}.jpg"
            have = any(e["year"] == year for e in entries)
            if frame_path.exists() and have and not args.force:
                print(f"{year}: exists, skipping", flush=True)
                continue
            t0 = time.time()
            try:
                result = process_year(catalog, year, aoi_poly, transform, width, height)
            except Exception as exc:
                print(f"{year}: ERROR {exc}", flush=True)
                continue
            if result is None:
                print(f"{year}: no usable scenes found", flush=True)
                continue
            Image.fromarray(result["img"]).save(frame_path, quality=88)
            entries = [e for e in entries if e["year"] != year] + [result["meta"]]
            write_manifest(out_dir, entries)
            m = result["meta"]
            print(
                f"{year}: {m['scenes']} scene(s) {m['dates']} {m['platforms']} "
                f"cloud={m['cloud']}% lake_cloud={m['aoi_cloud']} cov={m['coverage']} "
                f"area={m['area_km2']} km2 ({time.time() - t0:.0f}s)",
                flush=True,
            )
    print("done", flush=True)


if __name__ == "__main__":
    main()
