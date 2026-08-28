#%% Shock-front propagation streamlines
# Builds a shock arrival-time map T(r, phi) from dt_shock_track.csv, then uses
# the eikonal relation to get the speed AND direction of the front:
#
#     v_sh = 1/|grad T|      n = grad T/|grad T|      grad T = (dT/dr, (1/r) dT/dphi)
#
# The front moves from early arrival times to late ones, so +grad T is the
# propagation direction. Streamlines trace n, coloured by the four-velocity
# u = gamma*beta.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata

CSV     = '/scratch/aripoll/athena_out/outputs/dt_shock_track.csv'
R_STAR  = 1.0            # stellar radius, code units
SIGMA   = 3              # smoothing of the arrival map, in cells
PHI_LO, PHI_HI = 0, 90   # wedge to plot, degrees CCW from +x
RMIN, RMAX = 0.15, 1.5   # radial span to draw
NGRID   = 220            # Cartesian grid resolution for the interpolated field
DENSITY = 1.6            # streamplot line density
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
# streamplot needs an evenly spaced Cartesian grid, so interpolate the polar
# field -- defined only on shocked cells -- onto one, then blank every grid
# point that falls outside the [RMIN, RMAX] x [PHI_LO, PHI_HI] wedge.
R2, PHI2 = np.meshgrid(r_grid, phi, indexing='ij')
xs = (R2 * np.cos(PHI2)).ravel()
ys = (R2 * np.sin(PHI2)).ravel()
uxs = (n_r * np.cos(PHI2) - n_t * np.sin(PHI2)).ravel()   # polar -> Cartesian
uys = (n_r * np.sin(PHI2) + n_t * np.cos(PHI2)).ravel()
us  = u.ravel()

good = np.isfinite(uxs) & np.isfinite(uys)
pts  = np.column_stack([xs[good], ys[good]])

gx = np.linspace(0.0, RMAX, NGRID)
gy = np.linspace(0.0, RMAX, NGRID)
GX, GY = np.meshgrid(gx, gy)
UX = griddata(pts, uxs[good], (GX, GY), method='linear')
UY = griddata(pts, uys[good], (GX, GY), method='linear')
C  = griddata(pts, us[good],  (GX, GY), method='linear')

GR   = np.hypot(GX, GY)
GPHI = np.degrees(np.arctan2(GY, GX))
off  = (GR < RMIN) | (GR > RMAX) | (GPHI < PHI_LO) | (GPHI > PHI_HI)
UX[off], UY[off] = np.nan, np.nan

ok = np.isfinite(UX) & np.isfinite(C)
print(f'drawing streamlines over {ok.sum()} of {UX.size} grid points')

fig, ax = plt.subplots(figsize=(9, 7.5))
arc = np.linspace(PHI_LO, PHI_HI, 400) * np.pi / 180
ax.plot(R_STAR * np.cos(arc), R_STAR * np.sin(arc), color='grey', lw=1.4,
        label=r'$R = R_{*}$')

lo, hi = CLIM if CLIM else (np.nanmin(C[ok]), np.nanmax(C[ok]))
strm = ax.streamplot(gx, gy, UX, UY, color=C, cmap='plasma',
                     norm=LogNorm(vmin=max(lo, 1e-3), vmax=hi),
                     density=DENSITY, linewidth=1.2, arrowsize=1.0)
q = strm.lines

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
