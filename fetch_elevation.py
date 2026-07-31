"""Fetch Lake Powell daily pool elevation + storage from USBR hydrodata (site 919)
and write yearly values (mean of June-September, matching the imagery season)
to elevation.js / elevation.json."""

import json
import urllib.request
from collections import defaultdict
from pathlib import Path

BASE = "https://www.usbr.gov/uc/water/hydrodata/reservoir_data/919/csv"
DATATYPES = {"elevation_ft": 49, "storage_af": 17}


def fetch_daily(dt_id):
    url = f"{BASE}/{dt_id}.csv"
    with urllib.request.urlopen(url, timeout=60) as r:
        text = r.read().decode()
    out = {}
    for line in text.splitlines()[1:]:
        parts = line.strip().split(",")
        if len(parts) == 2 and parts[1]:
            out[parts[0]] = float(parts[1])
    return out


def summer_mean(daily, year):
    vals = [
        v
        for d, v in daily.items()
        if d.startswith(str(year)) and 6 <= int(d[5:7]) <= 9
    ]
    return round(sum(vals) / len(vals), 1) if vals else None


def main():
    out_dir = Path(__file__).parent
    series = {name: fetch_daily(dt) for name, dt in DATATYPES.items()}
    years = defaultdict(dict)
    for year in range(1963, 2027):
        for name, daily in series.items():
            v = summer_mean(daily, year)
            if v is not None:
                years[year][name] = v
    latest_date = max(series["elevation_ft"])
    data = {
        "latest": {
            "date": latest_date,
            "elevation_ft": series["elevation_ft"][latest_date],
        },
        "full_pool_ft": 3700.0,
        "years": {str(y): v for y, v in sorted(years.items())},
    }
    (out_dir / "elevation.json").write_text(json.dumps(data, indent=1))
    (out_dir / "elevation.js").write_text("const ELEVATION = " + json.dumps(data) + ";\n")
    print(f"wrote elevation for {len(years)} years, latest {latest_date}")


if __name__ == "__main__":
    main()
