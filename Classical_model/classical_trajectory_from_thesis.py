"""
File: classical_trajectories_1.py
Project: HHG-SaDAS
Code Description:
    | ~~~~~~~~~~~~: [Figure developed for my thesis] :~~~~~~~~~~~
    | Script for simulating classical electron trajectories in a
    | monochromatic laser field, relevant to the three-step model
    | of High Harmonic Generation (HHG). The code calculates and
    | visualizes:
    |   - Electron position and velocity as functions of time
    |   - Return times and recombination events
    |   - Kinetic energy at recombination relative to the
    |     ponderomotive potential (U_p)
    |
    | The simulation uses atomic units, laser parameters
    | (wavelength, intensity, frequency), and scans different
    | tunneling phases φ₀ to determine the maximum return energy.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Laser field: linearly polarized, continuous wave with finite cycles.
- Trajectories are obtained from classical equations of motion under
  the dipole approximation.
- Maximum return energy is expected to be ~3.17 U_p in agreement with
  the three-step model of HHG.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import numpy as np
from Atomic_units import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1 import make_axes_locatable
import Assistant.Decorate_axes.decorate_axes_L as da
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ---------------------- Common Figure Settings ----------------------
width = 6.2                                     # Width in inches
height = 6.0                                    # Height in inches
fig_scale_factor = 1.5                          # big=2 ; medium=1.5; small=1
tickslabel_size = 14
label_fontsize = 14
fig_size = (fig_scale_factor * width, fig_scale_factor * height)

dec_color = np.concatenate((da.mc.C_L, da.mc.des_col_1))
plt.rc('font', **{'family': 'serif', 'size': tickslabel_size})
da.mpl.rcParams['axes.prop_cycle'] = da.mpl.cycler(color=dec_color)
# -------------------------------------------------------------------

fig1 = plt.figure(figsize=fig_size)
ax1 = fig1.add_subplot(221); ax11 = fig1.add_subplot(222)
ax2 = fig1.add_subplot(223); ax22 = fig1.add_subplot(224)
ax2.set_axisbelow(True)             # <-- Ensure grid is below data
ax22.set_axisbelow(True)            # <-- Ensure grid is below data

axins1 = inset_axes(ax1, width="55%", height="30%", loc='lower left', borderpad=1.2)
axins1.set_yticklabels([])
axins1.tick_params(labelsize=10)

axins11 = inset_axes(ax11, width="55%", height="30%", loc='upper left', borderpad=1.2)
axins11.set_yticklabels([])
axins11.tick_params(labelsize=10)

da.decorate_2d([ax1, ax11, ax2, ax22])
da.decorate_2d([axins1, axins11], grid=False)


def E_field(t):
    return E0_au * np.cos(w0 * t)


def v(t, phi_0):
    return -(E0_au / w0) * (np.sin(w0 * t) - np.sin(phi_0))


def x(t, phi_0):
    return (E0_au / w0 ** 2) * (np.cos(phi_0) - np.cos(w0 * t) - (w0 * t - phi_0) * np.sin(phi_0))


def plot_x(t, phi_0):
    return x(t, phi_0) / (E0_au / w0 ** 2)


def KE(t, phi_0):
    return 0.5 * v(t, phi_0) ** 2


def sign_change_indices(arr):
    arr = arr[5:]                       # skip first few points
    return np.where(np.sign(arr[:-1]) != np.sign(arr[1:]))[0]



lambda_nm = 1064                               # wavelength (nm)
I0 = 5 * 10**13                                # Intensity (W/cm2)
I0_au = I0 / Int_0                             # Intensity (a.u)
E0_au = np.sqrt(I0_au)                         # Field intensity (a.u)
w0 = omega_au(lambda_nm); T = 2 * np.pi / w0
cpp = 2; tf = cpp * T
U_p = E0_au ** 2 / (4 * w0 ** 2)

phi_0_1st_qtr = np.linspace(0.0001, 0.99 * np.pi / 2, 60)               # first quarter period : φ0 = (0, π/4)
phi_0_3rd_qtr = np.linspace(1.001 * np.pi, 0.99 * 3 * np.pi / 2, 50)         # third quarter period : φ0 = (π, 3π/2)


def process_phase_range(phi_0_array):
    """Process a phase range and return results"""
    roots = []
    KE_max_return = []

    for i in range(len(phi_0_array)):
        t0 = phi_0_array[i] / w0
        t = np.linspace(t0, t0 + tf, 2000)

        x_t = x(t, phi_0_array[i])
        root_ind = sign_change_indices(x_t)[0]
        root = fsolve(lambda tt: x(tt, phi_0_array[i]), t[root_ind])[0]
        roots.append(root)
        KE_max_return.append(KE(root, phi_0_array[i]))

    KE_max_return = np.array(KE_max_return)
    return roots, KE_max_return


roots_low, KE_low = process_phase_range(phi_0_1st_qtr)
roots_high, KE_high = process_phase_range(phi_0_3rd_qtr)

# Combining all KE values for global normalization
all_KE = np.concatenate([KE_low, KE_high])
KE_min_global = all_KE.min()
KE_max_global = all_KE.max()

# Normalizing both ranges using global min/max
norm_low = (KE_low - KE_min_global) / (KE_max_global - KE_min_global + 1e-12)
norm_high = (KE_high - KE_min_global) / (KE_max_global - KE_min_global + 1e-12)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~: Plotting :~~~~~~~~~~~~~~~~~~~~~~~~~~
cmap = plt.get_cmap('jet')
print("Plotting first quarter period trajectories... φ0=(0, π/4)")
for i in range(len(phi_0_1st_qtr)):
    t0 = phi_0_1st_qtr[i] / w0
    t = np.linspace(t0, t0 + tf, 1500)
    x_t = plot_x(t, phi_0_1st_qtr[i])

    color = cmap(norm_low[i])
    ax1.plot(t / T, x_t, lw=1.5, color=color, alpha=0.8)
    axins1.plot(t / T, x(t, phi_0_1st_qtr[i]), lw=1.5, color=color, alpha=0.8)

print("Plotting third quarter period trajectories... φ0=(π, 3π/2)")
for i in range(len(phi_0_3rd_qtr)):
    t0 = phi_0_3rd_qtr[i] / w0
    t = np.linspace(t0, t0 + tf, 1500)
    x_t = plot_x(t, phi_0_3rd_qtr[i])

    color = cmap(norm_high[i])
    ax11.plot(t / T, x_t, lw=1.5, color=color, alpha=0.8)
    axins11.plot(t / T, x(t, phi_0_3rd_qtr[i]), lw=1.5, color=color, alpha=0.8)

for i in range(len(phi_0_3rd_qtr)):
    color = cmap(norm_high[i])
    axins11.scatter(roots_high[i] / T, x(roots_high[i], phi_0_3rd_qtr[i]), color=color, s=15)

KE_low_norm = KE_low / U_p
KE_high_norm = KE_high / U_p

KE_max_low_phi_0 = phi_0_1st_qtr[np.argmax(KE_low_norm)] / np.pi
KE_max_high_phi_0 = phi_0_3rd_qtr[np.argmax(KE_high_norm)] / np.pi

sc = ax2.scatter(phi_0_1st_qtr / np.pi, KE_low_norm,
                 c=KE_low_norm,        # color values based on data
                 cmap='jet',
                 label=rf'K$_{{max}}$ = {np.round(max(KE_low_norm), 3)}U$_p$',
                 zorder=2
                 )
ax2.axvline(float(KE_max_low_phi_0),
            lw=1.5,
            linestyle='dashed',
            color='gray',
            label=rf'φ$_0$={np.round(KE_max_low_phi_0, 3)}π',
            zorder=1
            )

ax22.scatter(phi_0_3rd_qtr / np.pi, KE_high_norm,
             c=KE_high_norm,        # color values based on data
             cmap='jet',
             label=rf'K$_{{max}}$ = {np.round(max(KE_high_norm), 3)}U$_p$',
             zorder=2
             )

ax22.axvline(float(KE_max_high_phi_0),
             lw=1.5,
             linestyle='dashed',
             color='gray',
             label=rf'φ$_0$={np.round(KE_max_high_phi_0, 3)}π',
             zorder=1
             )

divider = make_axes_locatable(ax22)
cax = divider.append_axes("right", size="3%", pad=0.1)  # adjust pad to shift right
cbar = fig1.colorbar(sc, cax=cax)
cbar.set_label(r'K$_{max}$ (a.u.)', fontsize=13)




imax_low = np.argmax(KE_low)
t0_max_low = phi_0_1st_qtr[imax_low] / w0
t_max_low = np.linspace(t0_max_low, t0_max_low + tf, 1500)
x_max_low = x(t_max_low, phi_0_1st_qtr[imax_low]) / (E0_au/w0**2)

ax1.annotate('', xy=(t_max_low[-50]/T, x_max_low[-50]), xytext=(t_max_low[-1]/T, x_max_low[-1]), arrowprops=dict(facecolor='blue', shrink=0.9, width=2, headwidth=4))


imax_high = np.argmax(KE_high)
t0_max_high = phi_0_3rd_qtr[imax_high] / w0
t_max_high = np.linspace(t0_max_high, t0_max_high + tf, 1500)
x_max_high = x(t_max_high, phi_0_3rd_qtr[imax_high]) / (E0_au/w0**2)

ax11.annotate('', xy=(t_max_high[-50]/T, x_max_high[-50]), xytext=(t_max_high[-1]/T, x_max_high[-1]), arrowprops=dict(facecolor='red', shrink=0.9, width=2, headwidth=4))

ax2.legend(loc='upper right', fontsize=12, framealpha=0.8, edgecolor='w')
ax22.legend(loc='upper right', fontsize=12, framealpha=0.8, edgecolor='w')

ax1.set_xlabel('t/T', fontsize=15)
ax11.set_xlabel('t/T', fontsize=15)
ax2.set_xlabel(r'φ$_0$/π', fontsize=label_fontsize)
ax22.set_xlabel(r'φ$_0$/π', fontsize=label_fontsize)
ax1.set_ylabel(r'z(t)/(E$_0$/w$_0^2$)', fontsize=15)
ax2.set_ylabel(r'KE/U$_p$', fontsize=label_fontsize)


sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=KE_min_global / U_p, vmax=KE_max_global / U_p))
divider = make_axes_locatable(ax11)
cax = divider.append_axes("right", size="3%", pad=0.1)
cbar = fig1.colorbar(sm, cax=cax)
cbar.set_label(r'K$_{max}$/U$_p$', fontsize=13)

axins1.axis([-0.005, 0.05, -0.005, 0.15])
axins11.axis([0.72, 0.91, -0.07, 0.06])


fig1.subplots_adjust(
    top=0.975,
    bottom=0.11,
    left=0.08,
    right=0.93,
    hspace=0.2,
    wspace=0.125
)

print(f"Global KE range: {KE_min_global / U_p:.3f} to {KE_max_global / U_p:.3f} U_p")

plt.show()