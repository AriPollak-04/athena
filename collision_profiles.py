"""Profiles of the antipodal collision outside the star, R = 1 to 1.2 R*.

Left column, profiles in phi at fixed R, centred on the collision axis. Right
column, profiles in R along the collision axis (solid) and off-axis (dashed).
One colour per snapshot. Run from the repo root.

    python collision_profiles.py jb_exp.out1.01500.athdf jb_exp.out1.01933.athdf

Velocities are four-velocity components (see CLAUDE.md). r0 is the jet tracer,
so r0 ~ 0 means stellar material.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "vis/python")
import athena_read  # noqa: E402

PHI_C = 90.0         # collision axis in degrees, +y in the slice plots
R_FIX = 1.1          # radius of the phi profiles
DPHI = 60.0          # half-width of the phi window, also the off-axis offset
R_LIM = (1.0, 1.2)
ROWS = [("rho", r"$\rho$", True), ("press", r"$P$", True),
        ("theta", r"$\Theta = P/\rho$", True),
        ("ur", r"$u^r = \gamma\beta_r$", False),
        ("uphi", r"$u^\phi = \gamma\beta_\phi$", False),
        ("r0", r"jet tracer $r_0$", False)]


def load(f):
    d = athena_read.athdf(f, quantities=["rho", "press", "vel1", "vel2", "r0"],
                          dtype=np.float64)
    q = {k: d[k][0] for k in ("rho", "press", "r0")}  # drop the nx3 = 1 axis, (nphi, nr)
    q["ur"], q["uphi"] = d["vel1"][0], d["vel2"][0]
    q["theta"] = q["press"] / q["rho"]                 # kT / (m c^2)
    return d["Time"], np.degrees(d["x2v"]), d["x1v"], q


def nearest(a, v):
    return int(np.argmin(np.abs(a - v)))


def plot(files, out="collision_profiles.pdf"):
    fig, ax = plt.subplots(len(ROWS), 2, figsize=(7, 10), sharex="col")
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(files)))
    for f, c in zip(files, colors):
        t, phi, r, q = load(f)
        i = nearest(r, R_FIX)
        win = np.abs(phi - PHI_C) <= DPHI
        on, off = nearest(phi, PHI_C), nearest(phi, PHI_C - DPHI / 2)
        rw = (r >= R_LIM[0]) & (r <= R_LIM[1])
        for row, (k, lbl, log) in enumerate(ROWS):
            ax[row, 0].plot(phi[win], q[k][win, i], color=c, label=f"$t={t:.2f}$")
            ax[row, 1].plot(r[rw], q[k][on, rw], color=c)
            ax[row, 1].plot(r[rw], q[k][off, rw], color=c, ls="--")
            ax[row, 0].set_ylabel(lbl)
            for a in ax[row]:
                a.set_yscale("log" if log else "linear")
    ax[0, 0].set_title(rf"$R = {R_FIX}\,R_*$")
    ax[0, 1].set_title(rf"solid $\phi={PHI_C:.0f}^\circ$, dashed "
                       rf"$\phi={PHI_C - DPHI / 2:.0f}^\circ$")
    ax[0, 0].legend(fontsize=7)
    ax[-1, 0].set_xlabel(r"$\phi$ [deg]")
    ax[-1, 1].set_xlabel(r"$R\ [R_*]$")
    fig.tight_layout()
    fig.savefig(out)
    return fig


if __name__ == "__main__":
    plot(sys.argv[1:])
