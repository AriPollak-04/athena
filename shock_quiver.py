"""Shock-front propagation quiver plot.

Builds a shock arrival-time map T(r, phi) from the tracker CSVs, then uses the
eikonal relation to recover both the speed and the direction of the front:

    v_sh = 1 / |grad T|            direction  n = grad T / |grad T|

with grad T = (dT/dr, (1/r) dT/dphi) in polar coordinates. The front moves from
early arrival times toward late ones, so +grad T is the propagation direction.

Three figure styles, all on the same eikonal field:
  'field'   - arrows throughout the quadrant on Cartesian axes, theta measured
              counter-clockwise from +x, coloured by the four-velocity gamma*beta
              (default).
  'surface' - the same layout, but arrows only where the front crosses R = R_*.
  'wedge'   - polar wedge, coloured by log10(v_sh / v_star).

Physical check: with the jet along phi = 0, arrival times rise toward the
equator, so leaving the star the front sweeps toward increasing theta and
n_phi > 0. tangential_check() reports this.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from scipy.ndimage import gaussian_filter

__all__ = ['arrival_map_from_dt_track', 'arrival_map_from_shock_front',
           'eikonal_field', 'four_velocity', 'surface_field',
           'tangential_check', 'plot_shock_quiver', 'plot_surface_quiver',
           'plot_field_quiver']


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
# 4. Surface-arrival quiver (Cartesian, theta measured CCW from +x)
# ----------------------------------------------------------------------------
def four_velocity(v):
    """u = gamma*beta from a 3-speed in units of c. NaN where v >= c."""
    v = np.asarray(v, dtype=float)
    out = np.full(v.shape, np.nan)
    ok = np.isfinite(v) & (v > 0) & (v < 1)
    out[ok] = v[ok] / np.sqrt(1.0 - v[ok] ** 2)
    return out


def surface_field(T, r_grid, phi_rad, r_star=1.0, sigma=3.0, phi_range=(0, 90)):
    """Front speed and direction where the shock crosses r = r_star.

    Returns (phi, v_sh, n_r, n_phi) for the angles inside phi_range, sampled at
    the radial index nearest r_star. The eikonal field is evaluated on the full
    periodic map first so the phi derivative wraps correctly.
    """
    v, n_r, n_t = eikonal_field(T, r_grid, phi_rad, sigma=sigma)
    i = int(np.argmin(np.abs(r_grid - r_star)))

    phi_deg = np.rad2deg(phi_rad)
    sel = (phi_deg >= phi_range[0]) & (phi_deg <= phi_range[1])
    if not sel.any():
        raise ValueError(f'no angles in phi_range={phi_range}')
    return phi_rad[sel], v[i, sel], n_r[i, sel], n_t[i, sel]


def tangential_check(phi, n_t, verbose=True):
    """Physical sanity check: leaving the star the front should sweep toward
    increasing theta, i.e. n_phi > 0. Returns the fraction of angles that do."""
    ok = np.isfinite(n_t)
    frac = float(np.mean(n_t[ok] > 0)) if ok.any() else np.nan
    if verbose:
        verdict = 'OK' if frac > 0.5 else 'FAILED - arrows sweep the wrong way'
        print(f'+theta sweep: {frac:.1%} of surface angles have n_phi > 0'
              f'   median n_phi {np.nanmedian(n_t):+.3f}   [{verdict}]')
    return frac


def _colour_values(v, color, v_star):
    """Scalar carried by the arrow colours, plus its colourbar label."""
    if color == 'u4':
        return four_velocity(v), r'$u^{\mu}$ 4-velocity'
    if color == 'logv':
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.log10(v / v_star), r'$\log\left(v_{sh}/v_{*}\right)$'
    raise ValueError("color must be 'u4' or 'logv'")


def _setup_cartesian_ax(ax, r_star, phi_range, xlim, ylim, title):
    """Equal-aspect x/y axes in R_*, with the grey R = R_* arc."""
    arc = np.linspace(phi_range[0], phi_range[1], 400) * np.pi / 180.0
    ax.plot(r_star * np.cos(arc), r_star * np.sin(arc), color='grey', lw=1.4,
            label=r'$R = R_{*}$', zorder=1)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x\ [R_{*}]$')
    ax.set_ylabel(r'$y\ [R_{*}]$')
    ax.legend(loc='upper right', fontsize=9)
    if title:
        ax.set_title(title)
    return ax


def _draw_quiver(ax, x, y, ux, uy, C, clim, cmap, arrow_scale, clabel, width):
    ok = np.isfinite(C) & np.isfinite(ux) & np.isfinite(uy)
    if not np.any(ok):
        raise ValueError('no finite arrows to draw; check r_star / sigma / clim')
    lo, hi = clim if clim is not None else (np.nanmin(C), np.nanmax(C))
    q = ax.quiver(x[ok], y[ok], ux[ok], uy[ok], C[ok],
                  cmap=cmap, norm=mcolors.Normalize(lo, hi),
                  angles='xy', scale_units='xy', scale=arrow_scale,
                  width=width, zorder=2)
    cb = plt.colorbar(q, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label(clabel, fontsize=11)
    return q


def plot_surface_quiver(T, r_grid, phi_rad, *, r_star=1.0, sigma=3.0,
                        phi_range=(0, 90), color='u4', v_star=1.0,
                        clim=None, cmap='plasma', arrow_scale=4.0,
                        uniform_length=False, xlim=(-0.05, 1.5),
                        ylim=(-0.05, 1.5), title=None, ax=None, check=True):
    """Shock velocity where the front crosses r = r_star, on Cartesian axes.

    theta runs counter-clockwise from +x, so phi_range=(0, 90) fills the first
    quadrant. Arrows are rooted on the arc and point along the front normal;
    colour is the four-velocity u = gamma*beta ('u4') or log10(v/v_star)
    ('logv'). Arrow length tracks the colour unless uniform_length.
    """
    phi, v, n_r, n_t = surface_field(T, r_grid, phi_rad, r_star=r_star,
                                     sigma=sigma, phi_range=phi_range)
    if check:
        tangential_check(phi, n_t)
    C, clabel = _colour_values(v, color, v_star)

    x, y = r_star * np.cos(phi), r_star * np.sin(phi)
    ux = n_r * np.cos(phi) - n_t * np.sin(phi)
    uy = n_r * np.sin(phi) + n_t * np.cos(phi)
    mag = np.ones_like(C) if uniform_length else C
    ux, uy = ux * mag, uy * mag

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 7.5))
    _setup_cartesian_ax(ax, r_star, phi_range, xlim, ylim, title)
    q = _draw_quiver(ax, x, y, ux, uy, C, clim, cmap, arrow_scale, clabel, 0.004)
    return ax, q


def plot_field_quiver(T, r_grid, phi_rad, *, r_star=1.0, sigma=3.0,
                      phi_range=(0, 90), color='u4', v_star=1.0, clim=None,
                      cmap='plasma', arrow_scale=18.0, uniform_length=True,
                      n_r_arrows=26, n_phi_arrows=26, rmin=None, rmax=None,
                      xlim=(-0.05, 1.5), ylim=(-0.05, 1.5), title=None,
                      ax=None, check=True):
    """Shock propagation across the whole quadrant, in the same Cartesian layout.

    Same axes, arc and colouring as plot_surface_quiver, but the arrows sample
    the full (r, phi) field rather than sitting on the r = r_star arc. Length is
    uniform by default so a dense field stays legible; the colour carries the
    magnitude.
    """
    v, n_r, n_t = eikonal_field(T, r_grid, phi_rad, sigma=sigma)

    keep = np.ones(r_grid.size, dtype=bool)
    if rmin is not None:
        keep &= r_grid >= rmin
    if rmax is not None:
        keep &= r_grid <= rmax
    r_idx = np.nonzero(keep)[0]
    if r_idx.size == 0:
        raise ValueError(f'no radii in [{rmin}, {rmax}]')
    ri = r_idx[np.linspace(0, r_idx.size - 1,
                           min(n_r_arrows, r_idx.size)).astype(int)]

    phi_deg_all = np.rad2deg(phi_rad)
    in_wedge = np.nonzero(
        (phi_deg_all >= phi_range[0]) & (phi_deg_all <= phi_range[1]))[0]
    if in_wedge.size == 0:
        raise ValueError(f'no angles in phi_range={phi_range}')
    pj = in_wedge[np.linspace(0, in_wedge.size - 1,
                              min(n_phi_arrows, in_wedge.size)).astype(int)]

    if check:
        i_s = int(np.argmin(np.abs(r_grid - r_star)))
        tangential_check(phi_rad[in_wedge], n_t[i_s, in_wedge])

    R, PHI = np.meshgrid(r_grid[ri], phi_rad[pj], indexing='ij')
    sub = np.ix_(ri, pj)
    V, NR, NT = v[sub], n_r[sub], n_t[sub]
    C, clabel = _colour_values(V, color, v_star)

    x, y = R * np.cos(PHI), R * np.sin(PHI)
    ux = NR * np.cos(PHI) - NT * np.sin(PHI)
    uy = NR * np.sin(PHI) + NT * np.cos(PHI)
    mag = np.ones_like(C) if uniform_length else C
    ux, uy = ux * mag, uy * mag

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 7.5))
    _setup_cartesian_ax(ax, r_star, phi_range, xlim, ylim, title)
    q = _draw_quiver(ax, x, y, ux, uy, C, clim, cmap, arrow_scale, clabel, 0.0035)
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

    # --- surface arrivals must sweep toward +theta (jet along phi=0)
    T, r_grid, phi_rad = arrival_map_from_dt_track('fake_rel.csv')
    phi_s, v_s, nr_s, nt_s = surface_field(T, r_grid, phi_rad, r_star=1.0)
    frac = tangential_check(phi_s, nt_s)
    u_s = four_velocity(v_s)
    print(f'surface    : u = gamma*beta spans '
          f'[{np.nanmin(u_s):.2f}, {np.nanmax(u_s):.2f}]  '
          f'(fixture: 1.29 equator -> 3.04 pole)')

    assert frac > 0.9, 'surface arrows do not sweep toward +theta'
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
    p.add_argument('--rmax', type=float, default=None)
    p.add_argument('--rticks', type=float, nargs='*', default=None)
    p.add_argument('--t-range', type=float, nargs=2, default=None, metavar=('T0', 'T1'))
    p.add_argument('--sigma', type=float, default=3.0,
                   help='arrival-map smoothing, cells')
    p.add_argument('--arrows', type=int, nargs=2, default=[22, 22],
                   metavar=('NR', 'NPHI'))
    p.add_argument('--clim', type=float, nargs=2, default=None,
                   metavar=('LO', 'HI'))
    p.add_argument('--cmap', default=None, help='matplotlib colormap name')
    p.add_argument('--style', choices=['field', 'surface', 'wedge'],
                   default='field',
                   help="'field': arrows throughout the quadrant, Cartesian "
                        "axes (default); 'surface': only on the R=R_* arc; "
                        "'wedge': polar")
    p.add_argument('--color', choices=['u4', 'logv'], default='u4',
                   help="'u4': four-velocity gamma*beta; 'logv': log10(v/v_star)")
    p.add_argument('--r-star', type=float, default=1.0, dest='r_star')
    p.add_argument('--arrow-scale', type=float, default=None, dest='arrow_scale',
                   help='smaller = longer arrows (data units per unit magnitude)')
    p.add_argument('--uniform-length', action='store_true', default=None,
                   dest='uniform_length',
                   help='draw every arrow the same length (default for --style field)')
    p.add_argument('--scale-by-magnitude', action='store_false',
                   dest='uniform_length',
                   help='arrow length tracks the colour (default for --style surface)')
    p.add_argument('--rmin', type=float, default=None)
    p.add_argument('--scale-by-speed', action='store_true')
    p.add_argument('--title', default=None)
    p.add_argument('--dpi', type=int, default=150)
    p.add_argument('--selftest', action='store_true')
    a = p.parse_args(argv)

    if a.selftest:
        return _selftest()
    if a.cmap is None:
        a.cmap = 'turbo' if a.style == 'wedge' else 'plasma'
    if a.uniform_length is None:
        a.uniform_length = (a.style == 'field')
    if a.arrow_scale is None:
        a.arrow_scale = 18.0 if a.style == 'field' else 4.0
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

    if a.style == 'field':
        _, ax = plt.subplots(figsize=(9, 7.5))
        plot_field_quiver(
            T, r_grid, phi_rad, r_star=a.r_star, sigma=a.sigma,
            phi_range=tuple(a.phi_range), color=a.color, v_star=a.v_star,
            clim=a.clim, cmap=a.cmap, arrow_scale=a.arrow_scale,
            uniform_length=a.uniform_length, rmin=a.rmin, rmax=a.rmax,
            n_r_arrows=a.arrows[0], n_phi_arrows=a.arrows[1],
            title=a.title, ax=ax)
    elif a.style == 'surface':
        _, ax = plt.subplots(figsize=(9, 7.5))
        plot_surface_quiver(
            T, r_grid, phi_rad, r_star=a.r_star, sigma=a.sigma,
            phi_range=tuple(a.phi_range), color=a.color, v_star=a.v_star,
            clim=a.clim, cmap=a.cmap, arrow_scale=a.arrow_scale,
            uniform_length=a.uniform_length, title=a.title, ax=ax)
    else:
        _, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(7, 6.5))
        plot_shock_quiver(
            T, r_grid, phi_rad, v_star=a.v_star, phi_range=tuple(a.phi_range),
            rmax=a.rmax if a.rmax else 1.5, rticks=a.rticks, sigma=a.sigma,
            n_r_arrows=a.arrows[0], n_phi_arrows=a.arrows[1],
            vmin=a.clim[0] if a.clim else -0.5,
            vmax=a.clim[1] if a.clim else 1.5,
            cmap=a.cmap, scale_by_speed=a.scale_by_speed, title=a.title, ax=ax)
    plt.tight_layout()
    plt.savefig(a.out, dpi=a.dpi, bbox_inches='tight')
    print(f'wrote {a.out}')


if __name__ == '__main__':
    main()
