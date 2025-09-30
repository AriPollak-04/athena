#%% Lets look at shock_breakout speeds
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
import astropy.constants as c

# Load the csv file
data = pd.read_csv('/Users/aripollak/work/shock_breakout.csv')
data.set_index('angle_deg', inplace=True)

#%%
#Quick plot
plt.figure(figsize=(10, 6))
plt.plot(data.index, data['breakout_time'], marker='o', label='Shock Speed', color = 'black')
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
