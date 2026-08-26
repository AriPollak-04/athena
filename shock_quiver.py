"""Shock-front propagation quiver plot.

Builds a shock arrival-time map T(r, phi) from the tracker CSVs, then uses the
eikonal relation to recover both the speed and the direction of the front:

    v_sh = 1 / |grad T|            direction  n = grad T / |grad T|

with grad T = (dT/dr, (1/r) dT/dphi) in polar coordinates. The front moves from
early arrival times toward late ones, so +grad T is the propagation direction.

Arrows show n, colored by log10(v_sh / v_star).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from scipy.ndimage import gaussian_filter

__all__ = ['arrival_map_from_dt_track', 'arrival_map_from_shock_front',
           'eikonal_field', 'plot_shock_quiver']


# ----------------------------------------------------------------------------
# 1. Arrival-time maps
# ----------------------------------------------------------------------------
def arrival_map_from_dt_track(csv, t_range=None):
    """T(r, phi) from dt_shock_track.csv  [time, r, phi_deg, dlogS_dt].

    Arrival time of a cell = the time its entropy jump d(logS)/dt peaked.
    Cells never crossing the tracker threshold stay NaN.
    """
    df = pd.read_csv(csv)
    if t_range is not None:
        df = df[df['time'].between(*t_range)]
    if df.empty:
        raise ValueError(f'no rows in {csv} for t_range={t_range}')

    peak = df.loc[df.groupby(['r', 'phi_deg'])['dlogS_dt'].idxmax()]
    r_grid = np.sort(peak['r'].unique())
    phi_deg = np.sort(peak['phi_deg'].unique())
    T = (peak.pivot(index='r', columns='phi_deg', values='time')
             .reindex(index=r_grid, columns=phi_deg).to_numpy(dtype=float))
    return T, r_grid, np.deg2rad(phi_deg)


def arrival_map_from_shock_front(csv, nr=256, t_range=None):
    """T(r, phi) from shock_front.csv.

    Columns: time, angle_deg, r_shock, vr, v_tang, mushroom.

    This tracker stores front *radii* per angle per time, so the map is built by
    binning: each (t, angle, r_shock) marks radius r_shock as reached at time t,
    and each cell keeps the earliest such time.
    """
    df = pd.read_csv(csv)
    if t_range is not None:
        df = df[df['time'].between(*t_range)]
    if df.empty:
        raise ValueError(f'no rows in {csv} for t_range={t_range}')

    phi_deg = np.sort(df['angle_deg'].unique())
    r_edges = np.linspace(df['r_shock'].min(), df['r_shock'].max(), nr + 1)
    r_grid = 0.5 * (r_edges[:-1] + r_edges[1:])

    ri = np.searchsorted(r_edges, df['r_shock'].to_numpy(), 'right') - 1
    ri = np.clip(ri, 0, nr - 1)
    pj = np.searchsorted(phi_deg, df['angle_deg'].to_numpy())

    T = np.full((nr, phi_deg.size), np.inf)
    np.minimum.at(T, (ri, pj), df['time'].to_numpy())   # earliest arrival per cell
    T[np.isinf(T)] = np.nan
    return T, r_grid, np.deg2rad(phi_deg)


# ----------------------------------------------------------------------------
# 2. Eikonal speed + direction
# ----------------------------------------------------------------------------
def eikonal_field(T, r_grid, phi_rad, sigma=3.0, min_weight=0.3, r_floor=0.05):
    """Return (v_sh, n_r, n_phi) on the (r, phi) grid.

    T is smoothed first with a normalized (NaN-aware) Gaussian so the gradient
    is not dominated by cell-to-cell jitter in the arrival times. phi is treated
    as periodic, so pass the FULL 0..2pi map here even when plotting a wedge.
    """
    mask = np.isfinite(T)
    if sigma > 0:
        weight = gaussian_filter(mask.astype(float), sigma, mode=['nearest', 'wrap'])
        num = gaussian_filter(np.where(mask, T, 0.0), sigma, mode=['nearest', 'wrap'])
        with np.errstate(invalid='ignore', divide='ignore'):
            Ts = num / weight
        Ts[weight < min_weight] = np.nan
    else:
        Ts = T.copy()

    dphi = float(phi_rad[1] - phi_rad[0])
    dTdr = np.gradient(Ts, r_grid, axis=0)
    dTdphi = (np.roll(Ts, -1, axis=1) - np.roll(Ts, 1, axis=1)) / (2 * dphi)

    g_r = dTdr
    g_t = dTdphi / np.maximum(r_grid, r_floor)[:, None]   # physical phi-component
    g_mag = np.hypot(g_r, g_t)

    with np.errstate(invalid='ignore', divide='ignore'):
        v = 1.0 / g_mag
        n_r, n_t = g_r / g_mag, g_t / g_mag

    bad = ~np.isfinite(Ts) | (g_mag == 0) | (r_grid[:, None] < r_floor)
    for a in (v, n_r, n_t):
        a[bad] = np.nan
    return v, n_r, n_t


# ----------------------------------------------------------------------------
# 3. Plot
# ----------------------------------------------------------------------------
def _pi_ticks(ax, phi_range):
    """Angle ticks labelled as fractions of pi, like the reference figure."""
    from fractions import Fraction
    lo, hi = phi_range
    step = 45 if (hi - lo) <= 180 else 90
    ticks = np.arange(lo, hi + 1e-9, step)
    labels = []
    for t in ticks:
        f = Fraction(int(round(t)), 180).limit_denominator(24)
        if f == 0:
            labels.append('0')
        elif f == 1:
            labels.append(r'$\pi$')
        elif f.numerator == 1:
            labels.append(rf'$\frac{{1}}{{{f.denominator}}}\pi$')
        else:
            labels.append(rf'$\frac{{{f.numerator}}}{{{f.denominator}}}\pi$')
    ax.set_thetagrids(ticks, labels)


def plot_shock_quiver(T, r_grid, phi_rad, *, v_star=1.0, phi_range=(0, 90),
                      rmax=None, n_r_arrows=22, n_phi_arrows=22, sigma=3.0,
                      vmin=-0.5, vmax=1.5, cmap='turbo', arrow_scale=26,
                      theta_zero='N', theta_direction=1, title=None, ax=None,
                      pi_labels=True, scale_by_speed=False, grid_alpha=0.0,
                      r_labels_on_far_edge=True, rticks=None):
    """Polar wedge of shock propagation arrows colored by log10(v_sh / v_star).

    phi_range is in degrees and selects the wedge to display; the eikonal field
    is always computed on the full periodic map first.

    scale_by_speed draws arrow length proportional to v_sh instead of unit
    length; the direction is identical either way.
    """
    v, n_r, n_t = eikonal_field(T, r_grid, phi_rad, sigma=sigma)

    # --- subsample to a legible arrow density, then restrict to the wedge
    ri = np.linspace(0, r_grid.size - 1, min(n_r_arrows, r_grid.size)).astype(int)
    phi_deg_all = np.rad2deg(phi_rad)
    in_wedge = np.nonzero(
        (phi_deg_all >= phi_range[0]) & (phi_deg_all <= phi_range[1]))[0]
    if in_wedge.size == 0:
        raise ValueError(f'no angles in phi_range={phi_range}')
    n_keep = min(n_phi_arrows, in_wedge.size)
    pj = in_wedge[np.linspace(0, in_wedge.size - 1, n_keep).astype(int)]

    R, PHI = np.meshgrid(r_grid[ri], phi_rad[pj], indexing='ij')
    sub = np.ix_(ri, pj)
    V, NR, NT = v[sub], n_r[sub], n_t[sub]

    if rmax is not None:
        V = np.where(R <= rmax, V, np.nan)

    # polar unit vectors -> Cartesian arrow components
    U = NR * np.cos(PHI) - NT * np.sin(PHI)
    W = NR * np.sin(PHI) + NT * np.cos(PHI)
    if scale_by_speed:
        U, W = U * V, W * V

    with np.errstate(invalid='ignore', divide='ignore'):
        C = np.log10(V / v_star)

    ok = np.isfinite(C) & np.isfinite(U) & np.isfinite(W)

    if ax is None:
        _, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(7, 6.5))
    ax.set_theta_zero_location(theta_zero)
    ax.set_theta_direction(theta_direction)
    ax.set_thetamin(phi_range[0])
    ax.set_thetamax(phi_range[1])
    ax.set_rlim(0, rmax if rmax is not None else r_grid.max())
    if rticks is not None:
        ax.set_rticks(rticks)

    q = ax.quiver(PHI[ok], R[ok], U[ok], W[ok], C[ok],
                  cmap=cmap, norm=mcolors.Normalize(vmin, vmax),
                  scale=arrow_scale, width=0.006, pivot='tail')

    if pi_labels:
        _pi_ticks(ax, phi_range)
    if r_labels_on_far_edge:
        # matplotlib pins r-labels to the theta-min spine on a wedge axes and
        # ignores set_rlabel_position, so draw them along theta-max by hand.
        far = np.deg2rad(phi_range[1])
        rspan = ax.get_rmax() - ax.get_rmin()
        ax.set_yticklabels([])
        for t in ax.get_yticks():
            if ax.get_rmin() <= t <= ax.get_rmax():
                ax.text(far, t - 0.035 * rspan, f'{t:g}', rotation=90,
                        ha='right', va='center', fontsize=10,
                        transform=ax.transData)
    ax.grid(alpha=grid_alpha, lw=0.5)
    if title:
        ax.set_title(title, pad=18)

    cb = plt.colorbar(q, ax=ax, fraction=0.045, pad=0.11)
    cb.set_label(r'$\log\left(v_{sh}/v_{*}\right)$', fontsize=13)
    return ax, q


# ----------------------------------------------------------------------------
def _selftest():
    """Validate against make_fake.py's analytic ground truth.

    Note the eikonal returns the speed NORMAL to the front, v^2/sqrt(v^2+v'^2),
    which equals the radial speed only where the front is locally isotropic.
    """
    band_lo, band_hi = 0.4, 1.2

    # --- isotropic control: normal speed must equal the radial speed exactly
    T, r_grid, phi_rad = arrival_map_from_dt_track('fake_iso.csv')
    v, n_r, n_t = eikonal_field(T, r_grid, phi_rad, sigma=3.0)
    band = (r_grid > band_lo) & (r_grid < band_hi)
    e_iso = abs(np.nanmedian(v[band]) - 0.35) / 0.35
    print(f'isotropic  : v_rec {np.nanmedian(v[band]):.4f}  vs 0.3500'
          f'   rel err {e_iso:.2%}')
    print(f'             n_r {np.nanmedian(n_r[band]):+.4f} (want +1)'
          f'   |n_phi| {np.nanmedian(np.abs(n_t[band])):.4f} (want 0)')

    # --- anisotropic: compare to the analytic eikonal solution
    T, r_grid, phi_rad = arrival_map_from_dt_track('fake_dt_shock_track.csv')
    v, n_r, n_t = eikonal_field(T, r_grid, phi_rad, sigma=3.0)
    band = (r_grid > band_lo) & (r_grid < band_hi)

    v_eik = np.load('truth_v_eik.npy')
    ev = np.abs(np.nanmedian(v[band], axis=0) - v_eik) / v_eik
    enr = np.abs(np.nanmedian(n_r[band], axis=0) - np.load('truth_n_r.npy'))
    ent = np.abs(np.nanmedian(n_t[band], axis=0) - np.load('truth_n_t.npy'))
    print(f'anisotropic: v_norm  median rel err {np.nanmedian(ev):.2%}'
          f'   max {np.nanmax(ev):.2%}')
    print(f'             n_r     median abs err {np.nanmedian(enr):.4f}'
          f'   max {np.nanmax(enr):.4f}')
    print(f'             n_phi   median abs err {np.nanmedian(ent):.4f}'
          f'   max {np.nanmax(ent):.4f}')

    assert e_iso < 0.02, 'isotropic speed recovery failed'
    assert np.nanmedian(v[band]) > 0, 'speed must be positive'
    assert np.nanmedian(ev) < 0.05, 'anisotropic speed recovery failed'
    assert np.nanmedian(enr) < 0.05 and np.nanmedian(ent) < 0.05, 'direction failed'
    print('OK')


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('csv', nargs='?', help='dt_shock_track.csv or shock_front.csv')
    p.add_argument('--source', default='auto',
                   choices=['dt_track', 'shock_front', 'auto'])
    p.add_argument('--out', default='shock_quiver.png')
    p.add_argument('--v-star', type=float, default=1.0,
                   help='velocity normalisation for log10(v_sh/v_star); '
                        'default 1.0 = code units, i.e. log10(v_sh/c)')
    p.add_argument('--phi-range', type=float, nargs=2, default=[0, 90],
                   metavar=('LO', 'HI'), help='wedge in degrees (default 0 90)')
    p.add_argument('--rmax', type=float, default=1.5)
    p.add_argument('--rticks', type=float, nargs='*', default=None)
    p.add_argument('--t-range', type=float, nargs=2, default=None, metavar=('T0', 'T1'))
    p.add_argument('--sigma', type=float, default=3.0,
                   help='arrival-map smoothing, cells')
    p.add_argument('--arrows', type=int, nargs=2, default=[22, 22],
                   metavar=('NR', 'NPHI'))
    p.add_argument('--clim', type=float, nargs=2, default=[-0.5, 1.5],
                   metavar=('LO', 'HI'))
    p.add_argument('--cmap', default='turbo',
                   help="'turbo' (default) or 'jet' to match the paper")
    p.add_argument('--scale-by-speed', action='store_true')
    p.add_argument('--title', default=None)
    p.add_argument('--dpi', type=int, default=150)
    p.add_argument('--selftest', action='store_true')
    a = p.parse_args(argv)

    if a.selftest:
        return _selftest()
    if not a.csv:
        p.error('csv is required (or pass --selftest)')

    src = a.source
    if src == 'auto':
        src = 'shock_front' if 'shock_front' in a.csv else 'dt_track'
    loader = (arrival_map_from_dt_track if src == 'dt_track'
              else arrival_map_from_shock_front)
    T, r_grid, phi_rad = loader(a.csv, t_range=a.t_range)
    print(f'{src}: arrival map {T.shape} '
          f'({np.isfinite(T).mean():.1%} of cells shocked), '
          f'r [{r_grid.min():.3f}, {r_grid.max():.3f}]')

    _, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(7, 6.5))
    plot_shock_quiver(T, r_grid, phi_rad, v_star=a.v_star, phi_range=tuple(a.phi_range),
                      rmax=a.rmax, rticks=a.rticks, sigma=a.sigma,
                      n_r_arrows=a.arrows[0], n_phi_arrows=a.arrows[1],
                      vmin=a.clim[0], vmax=a.clim[1], cmap=a.cmap,
                      scale_by_speed=a.scale_by_speed, title=a.title, ax=ax)
    plt.tight_layout()
    plt.savefig(a.out, dpi=a.dpi, bbox_inches='tight')
    print(f'wrote {a.out}')


if __name__ == '__main__':
    main()
