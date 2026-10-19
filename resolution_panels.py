"""Pressure slices for a resolution study, one PDF per resolution.

Each PDF is one column (two snapshots stacked). Every column has the same panel
size and vertical layout, so join_pdfs.py can place them side by side. Only the
first column carries y labels and timestamps, and only the last carries the colorbar.

    from resolution_panels import zlim_of, plot_column
    runs = [files2400, files1800, files1200, files900, files600]
    zlim = zlim_of([f for files in runs for f in files])
    for i, files in enumerate(runs):
        plot_column(files, f"res_{i}.pdf", zlim, first=i == 0, last=i == len(runs) - 1)
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import yt

FIELD = ("gas", "pressure")
FONT = "STIXGeneral"  # Times-like, matches AAS/MNRAS body text
FS = 7                # pt at printed size
AX = 1.18             # panel side in inches, five columns fit a 7.1 in two-column page
DPI = 600             # image pixels per printed inch, the PDF default of 100 looks blocky
# Margins in inches. EDGE is the inner margin, so joined columns sit 2*EDGE apart.
TOP, BOT, GAP, EDGE = 0.2, 0.36, 0.04, 0.03
LEFT = 0.42                               # y ticks and label, first column only
CB_PAD, CB_W, CB_LBL = 0.04, 0.07, 0.40   # colorbar, last column only

mpl.rcParams.update({
    "font.family": FONT,
    "mathtext.fontset": "stix",
    "font.size": FS,
    "pdf.fonttype": 42,  # TrueType, journals reject the default Type 3 fonts
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})


def zlim_of(files):
    """Color range over every file, so all columns share one scale."""
    ext = [yt.load(f).all_data().quantities.extrema(FIELD) for f in files]
    return min(e[0] for e in ext), max(e[1] for e in ext)


def plot_column(files, out, zlim, first=False, last=False):
    left = LEFT if first else EDGE
    right = CB_PAD + CB_W + CB_LBL if last else EDGE
    w, h = left + AX + right, BOT + 2 * AX + GAP + TOP
    rects = [[left / w, (BOT + AX + GAP) / h, AX / w, AX / h],
             [left / w, BOT / h, AX / w, AX / h]]
    cbar_rect = [(left + AX + CB_PAD) / w, BOT / h, CB_W / w, (2 * AX + GAP) / h]
    fig = plt.figure(figsize=(w, h))
    axes = [fig.add_axes(r) for r in rects]

    dss = [yt.load(f) for f in files]
    for i, (ds, ax) in enumerate(zip(dss, axes)):
        p = yt.SlicePlot(ds, "z", FIELD)
        p.set_buff_size(round(AX * DPI))  # one image pixel per printed pixel
        # yt's own default math font is Computer Modern, override it to match the text
        p.set_font({"family": FONT, "math_fontfamily": "stix", "size": FS})
        p.set_cmap(FIELD, "inferno")
        p.set_zlim(FIELD, *zlim)
        p.set_colorbar_label(FIELD, r"$P$")
        p.set_xlabel(r"$x\ [R_*]$")
        p.set_ylabel(r"$y\ [R_*]$")
        if first:
            # top center lies on the dark floor gas above the star in every panel
            p.annotate_timestamp(x_pos=0.5, y_pos=0.95, corner=None,
                                 time_format=r"$t = {time:.2f}\,R_*/c$",
                                 text_args={"color": "white", "size": FS})
        plot = p.plots[FIELD]
        plot.figure, plot.axes = fig, ax
        plot.cax = fig.add_axes(cbar_rect)
        p.render()
        # yt always draws a colorbar, keep only one (bottom panel, last column)
        if i == 0 or not last:
            plot.cax.remove()

    # render() resizes the figure and moves the axes to yt's own layout, undo that
    fig.set_size_inches(w, h)
    for ax, r in zip(axes, rects):
        ax.set_position(r)
    if last:
        plot.cax.set_position(cbar_rect)

    nx, ny = dss[0].domain_dimensions[:2]
    axes[0].set_title(rf"${nx} \times {ny}$", fontsize=FS + 1)
    axes[0].set_xlabel("")
    axes[0].tick_params(labelbottom=False)
    # after render, otherwise yt overwrites them
    for ax in axes:
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])
        if not first:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)

    # no bbox_inches="tight", it would give columns different sizes
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
