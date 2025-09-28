"""
~ Spectra free_H.py

• plots dipole moment d(t).
• plots dipole acceleration d_A(t).
• plots power spectrum dipole moment d(t).
• plots power spectrum dipole acceleration d_A(t).

"""
import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))   # Ensure project root (HHG-SaDAS) is in sys.path

from parameters_and_functions import dt, tf, w0, T
from Assistant.Decorate_axes import decorate_axes_D as da
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.fft as ft
import pandas as pd
import numpy as np
this_dir = Path(__file__).resolve().parent

dip_mom_file = 'Evo_steps=88036_He(1s)_m=0_SAE-M1__L=20_kmax=50_N=200_rmax=200_Lmap=80_dt=0.1.xlsx'
file_dip_mom = this_dir / 'Time_evolution_data' / f'{dip_mom_file}'
file_dip_mom = pd.read_excel(file_dip_mom, header=None, skiprows=1).to_numpy().T

t = file_dip_mom[0]; E_t = file_dip_mom[1]; dipole_moment = file_dip_mom[2]

dipole_mom_FFT = ft.fft(dipole_moment)
freq_w = ft.fftfreq(len(t), d=dt) * 2 * np.pi                   # Angular frequency: w=[-oo, oo]
dip_mom_power_spectra = np.abs(dipole_mom_FFT / tf)**2          # Harmonic power : dipole mom
positive_freq_mask = freq_w > 0
freq_pos = freq_w[positive_freq_mask]
harmonic_order = freq_pos / w0                                          # Harmonic Order (w/w0)
dip_mom_power_spectra = dip_mom_power_spectra[positive_freq_mask]       # P(w) : dipole mom


# ~~~~~~~~: Conversion efficiency :~~~~~~~~~~
harmonic_energy = np.trapz(dip_mom_power_spectra, freq_pos)
laser_energy = np.trapz(E_t**2, t)

eta = harmonic_energy / laser_energy
print(f"Conversion efficiency: {eta:.10f}")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Plotting :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
da.decorate_2d([ax1, ax2])

ax1.plot(t/T, dipole_moment, lw=2, label='d(t)')
ax2.plot(harmonic_order, np.log10(dip_mom_power_spectra), label=r'log$_{10}$(P($\omega$)) ~ d(t)')

ax1.set_xticks(np.arange(0, 61, 10))
ax2.set_xticks(np.arange(1, int(max(harmonic_order)), 2))
ax2.set_yticks(np.arange(-20, -2, 4))

ax1.set_ylim(1.5 * min(dipole_moment), 1.5 * max(dipole_moment)); ax2.set_xlim(-1, 102)

ax1.set_xlabel(r't/T :$\longrightarrow$', fontsize=15)
ax1.set_ylabel(r'd(t) :$\longrightarrow$', fontsize=15)
ax2.set_xlabel(r'ω/ω$_0$ :$\longrightarrow$', fontsize=15)
ax2.set_ylabel(r'log$_{10}$(P($\omega$)) :$\longrightarrow$', fontsize=15)

ax1.legend(loc='upper right', fontsize=12, framealpha=0.5, edgecolor='k')
ax2.legend(loc='upper right', fontsize=12, framealpha=0.5, edgecolor='k')

fig.suptitle(dip_mom_file, fontsize=12)
fig.subplots_adjust(top=0.94, bottom=0.069, right=0.97, left=0.048, hspace=0.14)


plt.show()