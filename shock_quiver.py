#%% Shock-front propagation quiver
# Builds a shock arrival-time map T(r, phi) from dt_shock_track.csv, then uses
# the eikonal relation to get the speed AND direction of the front:
#
#     v_sh = 1/|grad T|      n = grad T/|grad T|      grad T = (dT/dr, (1/r) dT/dphi)
#
# The front moves from early arrival times to late ones, so +grad T is the
# propagation direction. Arrows show n, coloured by the four-velocity u = gamma*beta.
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter

CSV     = '/scratch/aripoll/athena_out/outputs/dt_shock_track.csv'
R_STAR  = 1.0            # stellar radius, code units
SIGMA   = 3              # smoothing of the arrival map, in cells
PHI_LO, PHI_HI = 0, 90   # wedge to plot, degrees CCW from +x
RMIN, RMAX = 0.15, 1.5   # radial span for the arrows
NR, NPHI = 40, 40        # arrow density (blocks, not sample points)
MIN_WEIGHT = 0.15        # drop smoothed cells with less real data than this nearby
CLIM    = None           # e.g. (1.0, 30); None = autoscale (log)

#%% Arrival-time map: when did the entropy jump peak in each cell?
df = pd.read_csv(CSV)
peak = df.loc[df.groupby(['r', 'phi_deg'])['dlogS_dt'].idxmax()]

r_grid  = np.sort(peak['r'].unique())
phi_deg = np.sort(peak['phi_deg'].unique())
phi     = np.deg2rad(phi_deg)
T = (peak.pivot(index='r', columns='phi_deg', values='time')
         .reindex(index=r_grid, columns=phi_deg).to_numpy(dtype=float))

print(f'{T.shape[0]} radii x {T.shape[1]} angles, '
      f'{np.isfinite(T).mean():.1%} of cells shocked')

#%% Smooth, then differentiate
# Normalised (NaN-aware) Gaussian: smooth data and mask with the same kernel so
# unshocked cells don't drag the result toward zero. phi wraps, r does not.
mask   = np.isfinite(T)
weight = gaussian_filter(mask.astype(float), SIGMA, mode=['nearest', 'wrap'])
num    = gaussian_filter(np.where(mask, T, 0.0), SIGMA, mode=['nearest', 'wrap'])
with np.errstate(invalid='ignore', divide='ignore'):
    Ts = num / weight
Ts[weight < MIN_WEIGHT] = np.nan    # too little real data nearby to trust

dphi   = phi[1] - phi[0]
dTdr   = np.gradient(Ts, r_grid, axis=0)
dTdphi = (np.roll(Ts, -1, axis=1) - np.roll(Ts, 1, axis=1)) / (2 * dphi)

g_r = dTdr
g_t = dTdphi / r_grid[:, None]          # physical phi-component
g   = np.hypot(g_r, g_t)

with np.errstate(invalid='ignore', divide='ignore'):
    v_sh = 1.0 / g                      # speed normal to the front, in c
    n_r, n_t = g_r / g, g_t / g         # unit propagation direction
    u = np.where((v_sh > 0) & (v_sh < 1), v_sh / np.sqrt(1 - v_sh**2), np.nan)

# Physical check: with the jet along phi=0 the front arrives late at the equator,
# so on the way out it should sweep toward increasing theta, i.e. n_phi > 0.
# Mask BEFORE comparing -- NaN > 0 is False, and a bool array cannot hold NaN,
# so nanmean() on the comparison would silently count unshocked angles as failures.
i_surf   = np.argmin(np.abs(r_grid - R_STAR))
in_wedge = (phi_deg >= PHI_LO) & (phi_deg <= PHI_HI)
row      = n_t[i_surf, in_wedge]
ok_row   = np.isfinite(row)
sweep    = np.mean(row[ok_row] > 0) if ok_row.any() else np.nan
print(f'+theta sweep at R*: {sweep:.0%} of {ok_row.sum()} shocked angles '
      f'({ok_row.mean():.0%} of the arc)  '
      f'[{"OK" if sweep > 0.5 else "WRONG WAY"}]')

# Where the front runs nearly tangential, 1/|grad T| is a *pattern* speed and can
# exceed c. Those cells become NaN in u and are not drawn -- report how many.
fin = np.isfinite(v_sh)
print(f'u: median {np.nanmedian(u):.1f}, max {np.nanmax(u):.1f};  '
      f'{np.mean(v_sh[fin] >= 1):.2%} of shocked cells superluminal (not drawn), '
      f'{int(np.nansum(u > 10))} cells with u > 10')

#%% Plot
# Block-average onto a coarse grid. Subsampling by index would mostly land on
# unshocked cells when coverage is patchy; averaging keeps any block with data.
ri = np.nonzero((r_grid >= RMIN) & (r_grid <= RMAX))[0]
pj = np.nonzero(in_wedge)[0]
nr_b, nphi_b = min(NR, ri.size), min(NPHI, pj.size)
ri = ri[:(ri.size // nr_b) * nr_b]
pj = pj[:(pj.size // nphi_b) * nphi_b]

def blocks(a):
    # mean over each block, ignoring NaN; blocks with no data stay NaN
    b = a[np.ix_(ri, pj)].reshape(nr_b, ri.size // nr_b, nphi_b, pj.size // nphi_b)
    with warnings.catch_warnings():     # all-NaN blocks are expected
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(b, axis=(1, 3))

with np.errstate(invalid='ignore'):
    R   = blocks(np.broadcast_to(r_grid[:, None], n_r.shape))
    PHI = blocks(np.broadcast_to(phi[None, :], n_r.shape))
    C, NRc, NTc = blocks(u), blocks(n_r), blocks(n_t)

norm = np.hypot(NRc, NTc)               # averaging unit vectors shortens them
NRc, NTc = NRc / norm, NTc / norm

x, y = R * np.cos(PHI), R * np.sin(PHI)
ux = NRc * np.cos(PHI) - NTc * np.sin(PHI)     # polar -> Cartesian
uy = NRc * np.sin(PHI) + NTc * np.cos(PHI)
ok = np.isfinite(C) & np.isfinite(ux)
print(f'drawing {ok.sum()} of {C.size} blocks')

fig, ax = plt.subplots(figsize=(9, 7.5))
arc = np.linspace(PHI_LO, PHI_HI, 400) * np.pi / 180
ax.plot(R_STAR * np.cos(arc), R_STAR * np.sin(arc), color='grey', lw=1.4,
        label=r'$R = R_{*}$')

lo, hi = CLIM if CLIM else (np.nanmin(C[ok]), np.nanmax(C[ok]))
q = ax.quiver(x[ok], y[ok], ux[ok], uy[ok], C[ok], cmap='plasma',
              norm=LogNorm(vmin=max(lo, 1e-3), vmax=hi),
              angles='xy', scale_units='xy', scale=18, width=0.0035)

ax.set_xlim(-0.05, RMAX)
ax.set_ylim(-0.05, RMAX)
ax.set_aspect('equal')
ax.set_xlabel(r'$x\ [R_{*}]$')
ax.set_ylabel(r'$y\ [R_{*}]$')
ax.set_title(rf'Shock propagation ({PHI_LO}-{PHI_HI}$^\circ$)')
ax.legend(loc='upper right', fontsize=9)
fig.colorbar(q, ax=ax, fraction=0.045, pad=0.03).set_label(
    r'$u^{\mu}$ 4-velocity', fontsize=11)
plt.tight_layout()
plt.show()
# %%
