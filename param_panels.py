"""Pressure slices of jb_exph and the five jb_param runs at one output, in 2x3 panels.

Every run writes out1 at the same dt, so one output number is the same time in all six.
Pick it on the command line, or change FRAME below:

    python param_panels.py 1500              # writes param_01500.png
    python param_panels.py 1500 fig.pdf

or from a notebook:

    from param_panels import plot_frame
    plot_frame(1500)
"""
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import yt
from mpl_toolkits.axes_grid1 import AxesGrid

OUTPUTS = "/scratch/aripoll/athena_out/outputs"
FRAME = 1500          # default output number, the NNNNN in <problem_id>.out1.NNNNN.athdf
FIELD = ("gas", "pressure")
FONT = "STIXGeneral"  # Times-like, matches AAS/MNRAS body text
FS = 11
DPI = 300

# Top row splits a fixed L_peak between rho and Gamma, bottom row changes the timing
RUNS = [
    ("jb_exph",         r"baseline, $\Gamma = 31$, $t_0 = t_{\rm stop} = 1.2$"),
    ("jb_param_hirho",  r"high $\rho$, low $\Gamma = 15.5$"),
    ("jb_param_hiG",    r"high $\Gamma = 62$, low $\rho$"),
    ("jb_param_hitjet", r"long $t_{\rm jet}$, $t_0 = t_{\rm stop} = 2.4$"),
    ("jb_param_lotjet", r"short $t_{\rm jet}$, $t_0 = t_{\rm stop} = 0.6$"),
    ("jb_param_lot0",   r"low $t_0 = 0.8$, $t_{\rm stop} = 1.2$"),
]

mpl.rcParams.update({
    "font.family": FONT,
    "mathtext.fontset": "stix",
    "font.size": FS,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})


def plot_frame(frame=FRAME, out=None):
    dss = [yt.load(f"{OUTPUTS}/{pid}.out1.{frame:05d}.athdf") for pid, _ in RUNS]
    times = [float(ds.current_time) for ds in dss]
    # a run with a different out1 dt would put this frame at another time
    assert max(times) - min(times) < 1e-3, dict(zip((pid for pid, _ in RUNS), times))
    ext = [ds.all_data().quantities.extrema(FIELD) for ds in dss]
    zmin, zmax = min(e[0] for e in ext), max(e[1] for e in ext)

    fig = plt.figure(figsize=(10, 7))
    grid = AxesGrid(fig, 111, nrows_ncols=(2, 3), axes_pad=(0.08, 0.3), label_mode="L",
                    share_all=True, cbar_location="right", cbar_mode="single",
                    cbar_size="3%", cbar_pad=0.08)
    for i, ds in enumerate(dss):
        p = yt.SlicePlot(ds, "z", FIELD)
        p.set_buff_size(1000)  # yt's 800 px default looks blocky at DPI
        p.set_font({"family": FONT, "math_fontfamily": "stix", "size": FS})
        p.set_cmap(FIELD, "inferno")
        p.set_zlim(FIELD, zmin, zmax)
        p.set_colorbar_label(FIELD, r"$P$")
        p.set_xlabel(r"$x\ [R_*]$")
        p.set_ylabel(r"$y\ [R_*]$")
        if i == 0:  # the assert above makes it the time of every panel
            # the corner straddles the white outside r = 1.5 and the dark disk, box it
            p.annotate_timestamp(corner="upper_left", draw_inset_box=True,
                                 time_format=r"$t = {time:.2f}\,R_*/c$",
                                 text_args={"color": "white", "size": FS})
        plot = p.plots[FIELD]
        plot.figure, plot.axes, plot.cax = fig, grid[i].axes, grid.cbar_axes[i]
        p.render()

    # after render, otherwise yt overwrites them
    for ax, (_, title) in zip(grid, RUNS):
        ax.set_title(title, fontsize=FS)
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])

    fig.savefig(out or f"param_{frame:05d}.png", bbox_inches="tight", dpi=DPI)
    plt.show()


if __name__ == "__main__":
    plot_frame(int(sys.argv[1]) if len(sys.argv) > 1 else FRAME, *sys.argv[2:3])
