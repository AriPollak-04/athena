"""Join single-page PDFs side by side, left to right in the order given.

    python join_pdfs.py joined.pdf res_0.pdf res_1.pdf res_2.pdf res_3.pdf res_4.pdf

Needs pypdf (pip install pypdf). Text and lines stay vector.
"""
import sys

from pypdf import PdfReader, PdfWriter, Transformation


def join(out, paths):
    pages = [PdfReader(p).pages[0] for p in paths]
    heights = {round(float(p.mediabox.height), 1) for p in pages}
    if len(heights) != 1:
        raise ValueError(f"page heights differ {heights}, rows would not line up")
    writer = PdfWriter()
    canvas = writer.add_blank_page(width=sum(float(p.mediabox.width) for p in pages),
                                   height=heights.pop())
    x = 0.0
    for p in pages:
        box = p.mediabox
        canvas.merge_transformed_page(
            p, Transformation().translate(x - float(box.left), -float(box.bottom)))
        x += float(box.width)
    writer.write(out)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    join(sys.argv[1], sys.argv[2:])
