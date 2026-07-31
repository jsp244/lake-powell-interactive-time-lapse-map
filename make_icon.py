"""Generate lake.ico — generic lake icon: desert-tan tile, branching blue lake."""

from PIL import Image, ImageDraw

S = 256


def rounded_rect(draw, box, radius, fill):
    draw.rounded_rectangle(box, radius=radius, fill=fill)


def make_base():
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    # desert tile with subtle vertical shade
    rounded_rect(d, (8, 8, S - 8, S - 8), 52, (214, 168, 110, 255))
    shade = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    ds = ImageDraw.Draw(shade)
    rounded_rect(ds, (8, 8, S - 8, S - 8), 52, (255, 235, 200, 60))
    img = Image.alpha_composite(img, shade.crop((0, 0, S, S // 3)).resize((S, S)).point(lambda p: p))
    d = ImageDraw.Draw(img)

    # lake body: organic blob with canyon arms, drawn as overlapping shapes
    lake = (23, 84, 130, 255)
    lake_hi = (37, 106, 159, 255)
    # main body (lower left, like Wahweap/Padre Bay)
    d.ellipse((38, 128, 158, 216), fill=lake)
    d.ellipse((96, 116, 186, 194), fill=lake)
    # meandering channel to upper right
    pts = [
        (150, 150), (172, 132), (166, 112), (188, 100),
        (182, 82), (206, 72), (200, 52),
    ]
    for i in range(len(pts) - 1):
        x0, y0 = pts[i]
        x1, y1 = pts[i + 1]
        d.line((x0, y0, x1, y1), fill=lake, width=26)
        d.ellipse((x1 - 13, y1 - 13, x1 + 13, y1 + 13), fill=lake)
    # side arms
    d.line((120, 140, 92, 108), fill=lake, width=16)
    d.ellipse((84, 100, 100, 116), fill=lake)
    d.line((166, 112, 196, 118), fill=lake, width=14)
    d.ellipse((190, 111, 204, 125), fill=lake)
    d.line((70, 180, 46, 200), fill=lake, width=16)
    # subtle water highlight
    d.ellipse((60, 140, 130, 176), fill=lake_hi)
    d.ellipse((104, 130, 160, 168), fill=lake_hi)

    # wave accents on the main body
    for wy, wx0, wx1 in ((176, 66, 118), (192, 84, 140), (160, 118, 156)):
        d.arc((wx0, wy - 7, wx0 + 24, wy + 7), 200, 340, fill=(235, 244, 250, 230), width=5)
        d.arc((wx1 - 24, wy - 7, wx1, wy + 7), 200, 340, fill=(235, 244, 250, 230), width=5)
    return img


from pathlib import Path

out = Path(__file__).parent
img = make_base()
sizes = [16, 24, 32, 48, 64, 128, 256]
img.save(out / "lake.ico", sizes=[(s, s) for s in sizes])
img.resize((128, 128)).save(out / "lake_preview.png")
print("wrote lake.ico")
