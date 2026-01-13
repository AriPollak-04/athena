#%% Lets look at shock_breakout speeds
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
import astropy.constants as c



#%% lets look at the total energy and momentum conservation from the history file
hst_file = '/Users/aripollak/work/jet_blast.hst'
hst_data = pd.read_csv(hst_file, delim_whitespace=True, comment='#', header=None)
hst_data.columns = ['time', 'dt', 'mass', 'mom1', 'mom2', 'mom3',
                    'KE1', 'KE2', 'KE3', 'totE', 'Etot', 'Eexcess', 'Px_pos', 'Px_neg', 'Py_tot', 'Pz_tot', 'Pgas', 'Gamma_int', 'V_tot']


P_total = np.sqrt((hst_data['Px_pos'] + hst_data['Px_neg'])**2 + hst_data['Py_tot']**2 + hst_data['Pz_tot']**2)


v_eff = P_total / (hst_data['Etot'])
#%%

# Plot total energy conservation
plt.figure(figsize=(10, 6), dpi=300)

plt.scatter(hst_data['time'], hst_data['totE'], label='Total Energy (Etot)', 
          color='green', linestyle='--') 
#plt.plot(hst_data['time'], P_total, label='Total Momentum', color='red')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Total Energy Conservation')
plt.grid(True)
plt.legend()


plt.figure(figsize=(10, 6), dpi=300)
plt.scatter(hst_data['time'], (hst_data['Px_pos'] + hst_data['Px_neg']), label='P_x', color='purple', zorder=3)
plt.scatter(hst_data['time'], hst_data['Py_tot'], label='P_y', color='orange', zorder=2)
plt.scatter(hst_data['time'], hst_data['Pz_tot'], label='P_z', color='cyan', zorder=1)
plt.xlabel('Time')
plt.ylabel('Total Momentum')
plt.title('Total Momentum Conservation')
plt.grid(True)
plt.legend()

#%%

# Load the csv file
data = pd.read_csv('/Users/aripollak/work/shock_breakout.csv')
data.set_index('angle_deg', inplace=True)

#%%
#Quick plot
plt.figure(figsize=(10, 6))
plt.plot(data.index, data['breakout_time'], marker='o', color = 'black')
plt.xlabel('Angle (degrees)')
plt.ylabel('Breakout Time')
plt.title('Shock Breakout Time vs Angle')
plt.grid(True)

plt.legend()

# ct vs r*theta *******




data['r_theta'] = data['radius'] * np.deg2rad(data.index)
plt.figure(figsize=(10, 6))
plt.plot(data['breakout_time']*c.c, data['r_theta'], marker='o', label='Shock Speed', color = 'black')
plt.xlabel('Breakout Time (ct)')
plt.ylabel('r * theta (radians)')
plt.title('Shock Breakout Time vs r * theta')
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
import yt
from pathlib import Path
ds = yt.load('/Users/aripollak/work/jet_blast.out1.00079.athdf')
print(ds)
print(ds.current_time)
print(ds.field_list)        # native fields
print(ds.derived_field_list)
slc = yt.SlicePlot(ds, "z", ("gas", "pressure"))
slc.set_log(("gas","pressure"), True)
slc.show()
# %%
ts = yt.load('/Users/aripollak/work/jet_blast.out1.*.athdf')
frames_dir = Path.home() / "work" / "frames_density"
frames_dir.mkdir(parents=True, exist_ok=True)
print(len(ts))

for i, ds in enumerate(ts):
    p = yt.SlicePlot(ds, "z", ("gas", "density"))
    p.set_log(("gas", "density"), True)
    p.annotate_timestamp(corner="upper_left")
    # Optional: keep the color scaling fixed across time (avoids "breathing")
    # p.set_zlim(("gas","density"), 1e-6, 1e-2)
    p.save(frames_dir / f"frame_{i:05d}.png")
#%%
import subprocess

mp4_path = str(frames_dir / "density_z.mp4")
subprocess.run([
    "ffmpeg", "-y",
    "-framerate", "24",
    "-i", str(frames_dir / "frame_%05d.png"),
    "-c:v", "libx264",
    "-pix_fmt", "yuv420p",
    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
    mp4_path
], check=True)
from IPython.display import Video, display
display(Video(mp4_path, embed=True))

# %%
