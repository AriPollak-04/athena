#%% Lets look at shock_breakout speeds
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
import io

#%% lets look at the total energy and momentum conservation from the history file
hst_file = '/scratch/aripoll/athena_out/outputs/jet_blast.hst'
cols = ['time', 'dt', 'mass', 'mom1', 'mom2', 'mom3',
        'KE1', 'KE2', 'KE3', 'totE', 'scalar', 'Etot', 'Eexcess',
        'Px_pos', 'Px_neg', 'Py_tot', 'Pz_tot', 'Pgas', 'Gamma_int', 'V_tot']
with open(hst_file) as f:
    lines = f.readlines()
header_indices = [i for i, line in enumerate(lines) if line.startswith('# Athena++')]
runs = []
for i, start in enumerate(header_indices):
    end = header_indices[i+1] if i+1 < len(header_indices) else len(lines)
    df = pd.read_csv(io.StringIO(''.join(lines[start:end])), sep=r'\s+', comment='#', header=None, names=cols)
    runs.append(df)
hst_data = pd.concat(runs, ignore_index=True)


P_total = np.sqrt((hst_data['Px_pos'] + hst_data['Px_neg'])**2 + hst_data['Py_tot']**2 + hst_data['Pz_tot']**2)


v_eff = P_total / (hst_data['Etot'])
#%%

# Plot total energy conservation
plt.figure(figsize=(10, 6), dpi=300)
for i, run in enumerate(runs):
    plt.scatter(run['time'], run['totE'], label=f'Run {i+1}')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Total Energy Conservation')
plt.grid(True)
plt.legend()


markers = ['o', 's', '^']
plt.figure(figsize=(10, 6), dpi=300)
for i, run in enumerate(runs):
    m = markers[i % len(markers)]
    plt.scatter(run['time'], (run['Px_pos'] + run['Px_neg']), color='purple', marker=m, label=f'P_x Run {i+1}', zorder=3)
    plt.scatter(run['time'], run['Py_tot'], color='orange', marker=m, label=f'P_y Run {i+1}', zorder=2)
    plt.scatter(run['time'], run['Pz_tot'], color='cyan', marker=m, label=f'P_z Run {i+1}', zorder=1)
plt.xlabel('Time')
plt.ylabel('Total Momentum')
plt.title('Total Momentum Conservation')
plt.grid(True)
plt.legend()

#%%

# Load the csv file
data = pd.read_csv('/scratch/aripoll/athena_out/outputs/shock_breakout.csv')
data.set_index('angle_deg', inplace=True)

#%%
c = 1.0  # Speed of light in code units

#Quick plot
plt.figure(figsize=(10, 6))
plt.plot(data.index, data['breakout_time'], marker='o', color = 'black')
plt.xlabel('Angle (degrees)')
plt.ylabel('Breakout Time')
plt.title('Shock Breakout Time vs Angle')
plt.grid(True)

plt.legend()





plt.figure(figsize=(10, 6))
plt.plot(data['breakout_time']*c, data['radius'], marker='o', label='Shock Speed', color = 'black')
plt.xlabel('Breakout Time (ct)')
plt.ylabel(r'$R_{star}$')
plt.title('Shock Breakout Time vs Radius')
plt.grid(True)


# Now take the derivative of this curve to get the speed
# data['shock_speed'] = np.gradient(data['breakout_time'], data.index)
# plt.figure(figsize=(10, 6))
# plt.plot(data.index, data['shock_speed'], marker='o', label='Shock Speed', color = 'blue')
# plt.xlabel('Angle (degrees)')
# plt.ylabel('Shock Speed (dTime/dAngle)')
# plt.title('Shock Speed vs Angle')
# plt.grid(True)
# plt.legend()

# %%

#plot R vs angle with time as color in polar plot [0, 90]
plt.figure(figsize=(10, 10))
ax = plt.subplot(111, projection='polar')
sc = ax.scatter(np.deg2rad(data.index), data['radius'],
                c=data['breakout_time'], cmap='viridis', s=100)
ax.set_thetamin(0)
ax.set_thetamax(90)
plt.colorbar(sc, label='Breakout Time')
ax.set_xlabel('Radius')
plt.title('Shock Breakout Radius vs Angle')
plt.grid(True)
# %%
