#%% Shock-front propagation quiver
# Builds a shock arrival-time map T(r, phi) from dt_shock_track.csv, then uses
# the eikonal relation to get the speed AND direction of the front:
#
#     v_sh = 1/|grad T|      n = grad T/|grad T|      grad T = (dT/dr, (1/r) dT/dphi)
#
# The front moves from early arrival times to late ones, so +grad T is the
# propagation direction. Arrows show n, coloured by the four-velocity u = gamma*beta.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

CSV     = '/scratch/aripoll/athena_out/outputs/dt_shock_track.csv'
R_STAR  = 1.0            # stellar radius, code units
SIGMA   = 3              # smoothing of the arrival map, in cells
PHI_LO, PHI_HI = 0, 90   # wedge to plot, degrees CCW from +x
RMIN, RMAX = 0.15, 1.5   # radial span for the arrows
NR, NPHI = 26, 26        # arrow density
CLIM    = None           # e.g. (1.3, 3.0); None = autoscale to the data

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
Ts[weight < 0.3] = np.nan     # too little real data nearby to trust

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
i_surf = np.argmin(np.abs(r_grid - R_STAR))
in_wedge = (phi_deg >= PHI_LO) & (phi_deg <= PHI_HI)
sweep = np.nanmean(n_t[i_surf, in_wedge] > 0)
print(f'+theta sweep at R*: {sweep:.0%} of angles  '
      f'[{"OK" if sweep > 0.5 else "WRONG WAY"}]')

#%% Plot
ri = np.nonzero((r_grid >= RMIN) & (r_grid <= RMAX))[0]
ri = ri[np.linspace(0, ri.size - 1, min(NR, ri.size)).astype(int)]
pj = np.nonzero(in_wedge)[0]
pj = pj[np.linspace(0, pj.size - 1, min(NPHI, pj.size)).astype(int)]

R, PHI = np.meshgrid(r_grid[ri], phi[pj], indexing='ij')
sub = np.ix_(ri, pj)
C, NRc, NTc = u[sub], n_r[sub], n_t[sub]

x, y = R * np.cos(PHI), R * np.sin(PHI)
ux = NRc * np.cos(PHI) - NTc * np.sin(PHI)     # polar -> Cartesian
uy = NRc * np.sin(PHI) + NTc * np.cos(PHI)
ok = np.isfinite(C) & np.isfinite(ux)

fig, ax = plt.subplots(figsize=(9, 7.5))
arc = np.linspace(PHI_LO, PHI_HI, 400) * np.pi / 180
ax.plot(R_STAR * np.cos(arc), R_STAR * np.sin(arc), color='grey', lw=1.4,
        label=r'$R = R_{*}$')

lo, hi = CLIM if CLIM else (np.nanmin(C), np.nanmax(C))
q = ax.quiver(x[ok], y[ok], ux[ok], uy[ok], C[ok], cmap='plasma',
              clim=(lo, hi), angles='xy', scale_units='xy', scale=18,
              width=0.0035)

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
