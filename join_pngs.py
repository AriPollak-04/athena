"""Join PNGs side by side, left to right in the order given.

    python join_pngs.py -o joined.png a.png b.png c.png d.png

    d=/scratch/aripoll/athena_out/frames/pattern_speed
    python join_pngs.py -o $d/linked_joined.png \\
        $d/linked_00072.png $d/linked_00105.png $d/linked_00150.png $d/linked_00230.png

Needs Pillow (pip install pillow).
"""
import os
import sys

from PIL import Image


def join(out, paths):
    if os.path.abspath(out) in map(os.path.abspath, paths):
        raise ValueError(f"{out} is also an input, it would be overwritten")
    images = [Image.open(p).convert("RGBA") for p in paths]
    heights = {im.height for im in images}
    if len(heights) != 1:
        raise ValueError(f"image heights differ {heights}, rows would not line up")
    canvas = Image.new("RGBA", (sum(im.width for im in images), heights.pop()))
    x = 0
    for im in images:
        canvas.paste(im, (x, 0))
        x += im.width
    canvas.save(out)


if __name__ == "__main__":
    # -o is required so "join_pngs.py linked_*.png" can't overwrite linked_00000.png
    if len(sys.argv) < 4 or sys.argv[1] != "-o":
        sys.exit(__doc__)
    join(sys.argv[2], sys.argv[3:])
