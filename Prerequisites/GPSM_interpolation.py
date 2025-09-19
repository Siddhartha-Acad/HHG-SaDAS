"""
File: GPSM_interpolation.py
Project: HHG-SaDAS
Code Description:
    | *** [The generalised pseudospectral interpolation] ***
    | Demonstrates angular and radial interpolation. Results are visualized
    | through plots comparing exact values with their interpolated counterparts.
    | Detailed description is given in Appendix-G of the thesis.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import numpy as np
from pathlib import Path
from scipy.special import lpmv
import matplotlib.pyplot as plt
from scipy.special import legendre
import Assistant.Decorate_axes.decorate_axes_D as da
this_dir = Path(__file__).resolve().parent  # Relative file path system


def a_legendre(l, m, x):
    """
    Associated Legendre polynomial without Condon-Shortley phase
    """
    return lpmv(m, l, x) * (-1)**m

def P_n_prime(N, x):
    return (N*(N+1) / ((2*N+1)*(1-x**2))) * (legendre(N-1)(x) - legendre(N+1)(x))

def f(x):
    L_map = 80; r_max = 200
    alpha = 2 * L_map / r_max
    return L_map*(1 + x) / (1 - x + alpha)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Angular Interpolation :~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L = 10
l = 10; m = 0
theta_org = np.linspace(0, np.pi, 200)
roots, weights = np.polynomial.legendre.leggauss(L+1)
theta_k = np.arccos(roots); Ny = L+1

LegInterp = np.zeros(len(theta_org))
P_N_cos_theta = legendre(Ny)(np.cos(theta_org))


def ang_interpolation_funcn(theta):
    return a_legendre(l, m, np.cos(theta))

for j in range(Ny):
    LegInterp += ang_interpolation_funcn(theta_k[j]) * P_N_cos_theta / ((np.cos(theta_org) - np.cos(theta_k[j])) * P_n_prime(L+1, np.cos(theta_k[j])))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Radial Interpolation :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colloc_file = '199_AnaDeriv_collocation_points.txt'
colloc_pt_file = rf'E:\Python_programs\Pythonic_Physics\HHG\GPSM_Y_l0\Collocation_points\AnaDeriv_Colloc_pt\{colloc_file}'
colloc_pt = np.loadtxt(colloc_pt_file, skiprows=1, usecols=0)
x, dx = np.linspace(-1, 1, 500, retstep=True)
r = f(colloc_pt); Nx = len(colloc_pt)

g_tilde = (np.sin(4 * (np.pi/100) * r) + np.cos(3*(np.pi/100) * r)) * np.sin((np.pi / 200) * r)**2
g_Interp = np.zeros(len(x))

num_factor = (legendre(Nx+2)(x) - legendre(Nx)(x))
for i in range(Nx):
    g_Interp += g_tilde[i] * num_factor / ((2*Nx+3) * (x-colloc_pt[i]) * legendre(Nx+1)(colloc_pt[i]))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Plotting :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig1 = plt.figure()
fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(211)
ax3 = fig2.add_subplot(223)
ax4 = fig2.add_subplot(224)
da.decorate_2d([ax1, ax2, ax3, ax4])

ax1.plot(theta_k, ang_interpolation_funcn(theta_k), 'o', markersize=8, color='orange', label=rf'P$^{{{m}}}_{{{l}}}$(θ$_k$)', zorder=1)
ax1.plot(theta_org, ang_interpolation_funcn(theta_org), 'o-', lw=2, label=rf'P$^{{{m}}}_{{{l}}}$(θ) Exact', zorder=0)
ax1.plot(theta_org, LegInterp, '.-', lw=2, color='m', label=rf'P$^{{{m}}}_{{{l}}}$(θ) interpolated', zorder=0)
ax1.legend(loc='upper center', fontsize=15, framealpha=0.5, edgecolor='k')

ax2.plot(colloc_pt, g_tilde, 'o-', markersize=8, lw=2, label=r'$\tilde{g}$[r(x$_j$)]')
ax2.plot(x, g_Interp, 'o-', markersize=4, lw=2, color='m', label=r'$\tilde{g}$[r(x)] Interpolated')
ax2.legend(loc='lower left', fontsize=12, framealpha=0.5, edgecolor='k')

ax3.plot(np.diff(colloc_pt), 'o-', lw=2, label='dx (Gauss-Lobatto points)')
ax4.plot(np.diff(x), 'o-', lw=2, color='m', label='dx (Interpolating grid)'); ax4.set_ylim(-dx*1.5, dx*1.5)
ax3.legend(loc='lower center', fontsize=12, framealpha=0.5, edgecolor='k')
ax4.legend(loc='lower right', fontsize=12, framealpha=0.5, edgecolor='k')

plt.show()