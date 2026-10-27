"""Print-calibration strip for the printable's outline colour.

Why this exists
---------------
A printer reproduces a tint of ONE ink faithfully, because there is nothing to
mix.  A "light yellow" is not one ink: #ffd633 is yellow plus a little magenta,
so the printer has to mix two inks to build it, and when that mix is off the
sheet comes out pink or purple instead of yellow.  Every candidate here except
the flagged one is a tint of a single ink, so its hue has nothing to drift
towards.

(Blue is deliberately absent: cyan is a single ink but it reads as blue, and
cyan+magenta is a two-ink mix again.  Pure green would be cyan+yellow, two
inks, so it is out for the same reason as the yellow above.)

Print this at 100% / actual size and keep the row whose outline still looks
like itself.  Then pass it to the search:

    python3 qr_collider.py --emit out --verify --highlight '#ffb3ff'
"""

import os

from PIL import Image, ImageDraw, ImageFont

DPI = 300
SRC = os.path.join("out", "fillmap.png")
# the joined diagonal pair in the champion sheet, cropped for a live sample
CROP = (176, 208, 256, 288)

CANDIDATES = [
    ("#ffd9ff", "M 15%", "DEFAULT -- lightest magenta", False),
    ("#ffccff", "M 20%", "next step up if that looks faint", False),
    ("#ffb3ff", "M 30%", "clearly visible", False),
    ("#ff99ff", "M 40%", "strong", False),
    ("#ffff00", "Y 100%", "pure yellow -- needs yellow ink to be present", False),
    ("#ffff99", "Y 40%", "light yellow -- fades first on copies", False),
    ("#ffd633", "Y 80%+M 20%", "TWO INKS -- printed pink", True),
    ("#d9d9d9", "K 15%", "light grey -- neutral fallback", False),
]


def px(v):
    return max(1, int(round(v * DPI / 25.4)))


def font(size):
    for path in (r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\segoeui.ttf",
                 "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    try:
        return ImageFont.load_default(size)
    except TypeError:
        return ImageFont.load_default()


def sample(colour):
    """The real printable's outline, recoloured -- so the row shows the actual
    line weight and CORNER JOIN at print size, not a swatch in isolation."""
    if not os.path.exists(SRC):
        return None
    img = Image.open(SRC).convert("RGB").crop(CROP)
    d = img.load()
    for y in range(img.size[1]):
        for x in range(img.size[0]):
            r, g, b = d[x, y]
            if (r, g, b) != (255, 255, 255) and min(r, g, b) > 60:
                d[x, y] = colour          # ink pixel: not white, not module black
    return img.resize((img.size[0] * 2, img.size[1] * 2), Image.NEAREST)


def main():
    notes = ("A printer MIXES a light yellow from yellow + magenta, and that "
             "mix is what comes out pink.\nA tint of ONE ink has nothing to "
             "mix, so its hue cannot drift.\nEvery row below is such a tint, "
             "except the one flagged TWO INKS.\nPrint at 100% / actual size "
             "and keep the row whose outline still looks like itself.")
    row = px(15)
    pad = px(6)
    # start the rows below the whole header, whatever the header grows to
    top = px(12 + (notes.count("\n") + 1) * 4.9)
    width = px(150)
    height = top + row * len(CANDIDATES) + px(14)
    img = Image.new("RGB", (width, height), "white")
    dr = ImageDraw.Draw(img)
    f_big, f_small = font(px(5.4)), font(px(3.0))

    dr.text((pad, px(4)), "OUTLINE COLOUR CALIBRATION", fill=(0, 0, 0),
            font=f_big)
    dr.multiline_text((pad, px(11)), notes, fill=(0, 0, 0), font=f_small,
                      spacing=px(1.6))

    for i, (hexv, ink, note, flagged) in enumerate(CANDIDATES):
        y = top + i * row
        colour = tuple(int(hexv[1 + 2 * k:3 + 2 * k], 16) for k in range(3))
        dr.rectangle([pad, y, width - pad, y + row - px(2)],
                     outline=(200, 200, 200))
        dr.rectangle([pad + px(1), y + px(1), pad + px(20), y + row - px(3)],
                     fill=colour)
        label = f"{hexv}   {ink:10s}   {note}"
        dr.text((pad + px(23), y + px(1.6)), label, fill=(0, 0, 0), font=f_small)
        if flagged:
            dr.text((pad + px(23), y + px(6.6)), "^^ two inks mixed -- expect "
                    "pink or purple", fill=(0, 0, 0), font=f_small)
        s = sample(colour)
        if s is not None:
            img.paste(s, (width - pad - s.size[0], y + px(1.4)))

    dr.text((pad, height - px(9)), "sample at right = the real fill outline "
            "(the joined diagonal pair) in that colour, shown 2x", fill=(0, 0, 0),
            font=f_small)
    img.save("printer-swatches.png", dpi=(DPI, DPI))
    print(f"wrote printer-swatches.png  {img.size[0]}x{img.size[1]} px = "
          f"{img.size[0] / DPI * 25.4:.0f}x{img.size[1] / DPI * 25.4:.0f} mm "
          f"at {DPI} dpi")
    for hexv, ink, note, flagged in CANDIDATES:
        print(f"  {hexv}  {ink:14s} {note}")


if __name__ == "__main__":
    main()
