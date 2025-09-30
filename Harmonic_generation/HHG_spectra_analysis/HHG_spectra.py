"""
~ Spectra free_H.py

• plots dipole moment d(t).
• plots dipole acceleration d_A(t).
• plots power spectrum dipole moment d(t).
• plots power spectrum dipole acceleration d_A(t).

"""
import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import numpy as np
import pandas as pd
import scipy.fft as ft
from pathlib import Path
import matplotlib.pyplot as plt
from Assistant.Decorate_axes import decorate_axes_L_thesis as da
from Harmonic_generation.parameters_and_functions import dt, tf, w0, T

# ~~~~~~~~~~~~~~~~~~~~~~: Common Figure Settings :~~~~~~~~~~~~~~~~~~~~~
width = 6.2                         # Width in inches
height = 4                          # Height in inches
fig_scale_factor = 2                # big=2 ; medium=1.5; small=1
fig_size = (fig_scale_factor*width, fig_scale_factor*height)

plt.rc('font', **{'family': 'serif', 'size': 14})
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                Import and read time evolution data                 |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
this_dir = Path(__file__).resolve().parent

evo_data_file_name = 'Evo_steps=88036_He(1s)_m=0_SAE-M1__L=20_kmax=50_N=200_rmax=200_Lmap=80_dt=0.1.xlsx'
evo_data_file_path = this_dir.parent / 'Time_evolution_data' / f'{evo_data_file_name}'
evo_data = pd.read_excel(evo_data_file_path, header=None, skiprows=1).to_numpy().T

t = evo_data[0]                  # time               : (a.u.)
E_t = evo_data[1]                # Electric field     : (a.u.)
dipole_moment = evo_data[2]      # dipole moment d(t) : (a.u.)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    Calculating the HHG Spectra                     |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dipole_mom_FFT = ft.fft(dipole_moment)
freq_w = ft.fftfreq(len(t), d=dt) * 2 * np.pi                   # Nyquist Angular frequency (omega -> w).
dip_mom_power_spectra = np.abs(dipole_mom_FFT / tf)**2          # Harmonic power spectra.

positive_freq_mask = freq_w > 0             # Select positive angular frequencies; returns a boolean array where True corresponds to freq_w > 0
freq_pos = freq_w[positive_freq_mask]       # Positive Nyquist frequencies.

harmonic_order = freq_pos / w0                                          # Harmonic Order (w/w0)
dip_mom_power_spectra = dip_mom_power_spectra[positive_freq_mask]       # P(w) : The power spectra you want!


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Conversion efficiency                        |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
harmonic_energy = np.trapezoid(dip_mom_power_spectra, freq_pos)
laser_energy = np.trapezoid(E_t**2, t)

eta = harmonic_energy / laser_energy
print(f"Conversion efficiency: {eta:.10f}")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Plotting                              |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig = plt.figure(figsize=fig_size)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
da.decorate_2d([ax1, ax2])

ax1.plot(t/T, dipole_moment, color='#7F3FBF', lw=2, label='d(t)')
ax2.plot(harmonic_order, np.log10(dip_mom_power_spectra), color='#2CA02C', label=r'log$_{10}$(P($\omega$))')

ax1.set_xticks(np.arange(0, 61, 5))                                # Some controls
ax2.set_yticks(np.arange(-20, -2, 4))                               # over
ax2.set_xticks(np.arange(1, int(max(harmonic_order)), 4))           # dipole moment
ax2.set_xlim(-1, 102)                                               # and
ax1.set_ylim(1.5 * min(dipole_moment), 1.5 * max(dipole_moment))    # HHG spectra plot.

ax1.set_xlabel('Optical Cycles (t/T)', fontsize=15)
ax1.set_ylabel('d(t)', fontsize=15)
ax2.set_xlabel(r'Harmonic Order (ω/ω$_0$)', fontsize=15)
ax2.set_ylabel(r'log$_{10}$[P($\omega$)]', fontsize=15)

ax1.legend(loc='upper right', fontsize=14, framealpha=0.5, edgecolor='w')
ax2.legend(loc='upper right', fontsize=14, framealpha=0.5, edgecolor='w')

fig.suptitle(evo_data_file_name, fontsize=12)
fig.subplots_adjust(top=0.94, bottom=0.069, right=0.98, left=0.079, hspace=0.195)


plt.show()