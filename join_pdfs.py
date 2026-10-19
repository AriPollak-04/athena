"""Join single-page PDFs side by side, left to right in the order given.

    python join_pdfs.py -o joined.pdf res_0.pdf res_1.pdf res_2.pdf res_3.pdf res_4.pdf

Needs pypdf (pip install pypdf). Text and lines stay vector.
"""
import os
import sys

from pypdf import PdfReader, PdfWriter, Transformation


def join(out, paths):
    if os.path.abspath(out) in map(os.path.abspath, paths):
        raise ValueError(f"{out} is also an input, it would be overwritten")
    # cropbox, not mediabox, since pypdf clips merged pages to it (same if never cropped)
    pages = [PdfReader(p).pages[0] for p in paths]
    heights = {round(float(p.cropbox.height), 1) for p in pages}
    if len(heights) != 1:
        raise ValueError(f"page heights differ {heights}, rows would not line up")
    writer = PdfWriter()
    canvas = writer.add_blank_page(width=sum(float(p.cropbox.width) for p in pages),
                                   height=heights.pop())
    x = 0.0
    for p in pages:
        box = p.cropbox
        canvas.merge_transformed_page(
            p, Transformation().translate(x - float(box.left), -float(box.bottom)))
        x += float(box.width)
    writer.write(out)


if __name__ == "__main__":
    # -o is required so "join_pdfs.py res_*.pdf" can't overwrite res_0.pdf
    if len(sys.argv) < 4 or sys.argv[1] != "-o":
        sys.exit(__doc__)
    join(sys.argv[2], sys.argv[3:])
