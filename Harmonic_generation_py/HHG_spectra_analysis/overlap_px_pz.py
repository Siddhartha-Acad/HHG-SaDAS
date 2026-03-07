import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fig = plt.figure()
ax1 = fig.add_subplot(111)

file_px = 'VEvo_nopt=88036_Ne(2p)_m=1_SAE-M2_L=20_kmax=50_N=200_rmax=200_Lmap=80_dt=0.1.dat'
file_pz = 'VEvo_nopt=88036_Ne(2p)_m=0_SAE-M2_L=20_kmax=50_N=200_rmax=200_Lmap=80_dt=0.1.dat'

path = '/home/siddhartha/Documents/Academia/Programming/Python_programs/HHG-SaDAS/Harmonic_generation_py/Time_evolution_data'

file_pz = f'{path}/{file_pz}'
file_px = f'{path}/{file_px}'

# Read .dat files (assuming whitespace separated)
data_file_pz = pd.read_csv(file_pz, delim_whitespace=True, header=None, skiprows=1).to_numpy().T
data_file_px = pd.read_csv(file_px, delim_whitespace=True, header=None, skiprows=1).to_numpy().T


t = data_file_pz[0]
dipole_pz = data_file_pz[2]
dipole_px = data_file_px[2]

ax1.plot(t, dipole_pz, lw=2.5, label='d(t) ~ pz')
ax1.plot(t, dipole_px, lw=2.5, label='d(t) ~ px')

plt.legend()
plt.xlabel('Time (a.u.)')
plt.ylabel('Dipole moment (a.u.)')
plt.tight_layout()
plt.show()
