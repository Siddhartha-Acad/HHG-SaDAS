"""
File: classical_trajectories_1.py
Project: HHG-SaDAS
Code Description:
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
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Atomic_units import Int_0, omega_au
import Assistant.Decorate_axes.decorate_axes_L as da

fig1 = plt.figure(); fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
da.decorate_2d([ax1, ax2])

def E_field(t):
    return E0_au * np.cos(w0*t)
    # return E0_au * np.sin(w0*t)

def v(t, phi_0):
    return (E0_au/w0) * (np.sin(w0*t) - np.sin(phi_0))
    # return (E0_au/w0) * (np.cos(w0*t) - np.cos(phi_0))

def x(t, phi_0):
    return (E0_au/w0**2) * (np.cos(phi_0) - np.cos(w0*t) - (w0*t-phi_0) * np.sin(phi_0))
    # return (E0_au/w0**2) * (-np.sin(phi_0) + np.sin(w0*t) - (w0*t-phi_0) * np.cos(phi_0))

def KE(t, phi_0):
    return 0.5 * v(t, phi_0)**2

def sign_change_indices(arr):
    arr = arr[5:]
    return np.where(np.sign(arr[:-1]) != np.sign(arr[1:]))[0]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~: LASER info :~~~~~~~~~~~~~~~~~~~~~~~~~~
lambda_nm = 1064                            # wavelength (nm)
I0 = 5 * 10**13                             # Intensity (W/cm2)
I0_au = I0 / Int_0                          # Intensity (a.u)
E0_au = np.sqrt(I0_au)                      # Field intensity (a.u)
w0 = omega_au(lambda_nm); T = 2*np.pi / w0
cpp = 2; tf = cpp*T; dt = 1                 # cpp = cycles per pulse

# phi_0 = np.linspace(0.0001, 0.99*np.pi/2, 30)
# phi_0 = np.linspace(np.pi/2, np.pi, 30)                   # never recombines
phi_0 = np.linspace(1.001*np.pi, 0.99*3*np.pi/2, 50)
# phi_0 = np.linspace(3*np.pi/2,  2*np.pi, 50)              # never recombines


U_p = E0_au**2 / (4*w0**2)                                  # Ponderomotive force

KE_max_return = []
cmap = plt.get_cmap('jet')
colors = cmap(np.linspace(0, 1, len(phi_0)))
for i in range(len(phi_0)):
    t0 = phi_0[i] / w0
    t = np.linspace(t0, t0+tf, 2000)

    x_t = x(t, phi_0[i])
    root_ind = sign_change_indices(x_t)[0]
    root = fsolve(x, t[root_ind], phi_0[i])[0]
    KE_max_return.append(KE(root, phi_0[i]))

    if i == 0:
        label = rf'φ$_0$={np.round(phi_0[i]/np.pi, 3)}π'
    elif i == len(phi_0) - 1:
        label = rf'φ$_0$={np.round(phi_0[i]/np.pi, 3)}π'
    else:
        label = None

    ax1.plot(t/T, x_t, lw=1.5, color=colors[i], label=label)
    ax1.scatter(root/T, x(root, phi_0[i]), color=colors[i])

KE_max_return = np.array(KE_max_return) / U_p
KE_max_return_phi_0 = phi_0[np.where(KE_max_return==max(KE_max_return))[0]][0]/np.pi

ax2.axvline(KE_max_return_phi_0, lw=2, linestyle='dashed', color='royalblue', label=rf'φ$_0$={np.round(KE_max_return_phi_0, 3)}π')
ax2.plot(phi_0/np.pi, KE_max_return, 'o-', label=rf'max(KE) = {np.round(max(KE_max_return), 3)}U$_p$')
if any(label is not None for label in [1, len(phi_0) - 1]):
    ax1.legend(fontsize=15, framealpha=0.5, edgecolor='w')
ax2.legend(fontsize=15, framealpha=0.5, edgecolor='w')

ax1.set_title(r'Electron trajectory, varying tunneling time (t$_0$ = φ$_0$/ω$_0$)', fontsize=15)
ax2.set_title(r'Kinetic Energy at the first recombination', fontsize=15)
ax1.set_xlabel('t/T', fontsize=15); ax2.set_xlabel(r'φ$_0$/π', fontsize=15)
ax1.set_ylabel('z(t)', fontsize=15); ax2.set_ylabel(r'KE/U$_p$', fontsize=15)

plt.show()
