"""
File: Partial_wave_Expansion_m=0.py
Project: HHG-SaDAS
Code Description:
    | *** [VERY IMPORTANT CODE] ***
    | This script demonstrates the partial wave expansion of an arbitrary wavefunction particularly with m=0
    |
    | We consider the wavefunction:
    |     ψ(r, θ) = R(r) * Y_l0(cos θ)
    |
    | - The full ψ(r, θ) includes both radial and angular contributions.
    | - The corresponding partial waves only represent the radial function R(r).
    |
    | Example:
    | ---------
    | def Psi(r, theta, l=1):
    |     return np.sqrt((2*l + 1) / (4*np.pi)) * np.exp(-r**2) * legendre(l)(np.cos(theta))
    |
    | - Here, the true radial function is only exp(-r²).
    | - The pre-factor √((2l+1)/(4π)) and the Legendre polynomial are purely angular
    |   contributions.
    |
    | Key Insight:
    | ------------
    | - If we plot Psi(r, 0), it looks like a "radial distribution" at θ=0,
    |   but it actually contains l-dependent normalization factors.
    | - The pure radial function should only depend on r, without any l-dependent term.
    | - For example:
    |     If R(r) = 2 * exp(-r²), then g_l(r) will also contain the factor 2,
    |     but it will not include the l-dependent normalization factor.
    |
    | Notes:
    | ------
    | - This formulation is based on the 3D expansion of ψ(r, θ, φ).
    | - It specifically assumes m = 0, meaning the angular dependence is in terms of
    |   Legendre polynomials P_l(cos θ) rather than spherical harmonics Y_lm.
    | - For example:
    |     f(r, θ) = norm_const * exp(-r²) * P_l(cos θ)
    |   corresponds to the angular shape of the pz orbital (m = 0, dumbbell along z-axis).
    |
    | - However:
    | ----------
    |     f(r, θ) = norm_const * exp(-r²) * P_l(cos(θ + π/2))
    |   mimics the px orbital (m = 1, dumbbell along x-axis), but expanding it in terms
    |   of P_l(cos θ) will give the wrong function.
    |   Therefore, to approximate orbitals with m ≠ 0, a more general spherical harmonic
    |   expansion is required.
    |
    | More details in :
    | -----------------
    | Section-3.3.3 : 'Generalization over magnetic quantum number m'
    | Appendix-D: 'Derivation of the general partial wave formula'


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- This code is the first effort for calculating the partial waves gl(r) and their valudation.
- Don't forget, here the magnetic quantum number m = 0.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
from Assistant.Decorate_axes import decorate_axes_D as da

fig1 = plt.figure()
fig2 = plt.figure()
ax1 = fig1.add_subplot(111)             # g_l(r)
ax2 = fig2.add_subplot(121)             # psi(r, θ)
ax3 = fig2.add_subplot(122)             # psi(r, θ) expanded in Legendre Polynomial
da.decorate_2d([ax1])
da.decorate_polar([ax2, ax3])

def f(x):
    L_map = 1
    r_max = 5
    alpha = 2 * L_map / r_max
    return L_map*(1 + x) / (1 - x + alpha)

def Y_l0(l, theta):     # JUST try rotating anti-clk by adding extra pi/2 with theta, (i.e., theta -> theta+np.pi/2)
    return np.sqrt((2*l + 1) / (4*np.pi)) * legendre(l)(np.cos(theta))

def Psi(l, r, theta):
    return 2 * r * np.exp(-r**2) * Y_l0(l, theta)


l = 5
L = 5           # Number of angular collocation points: L+1   [NOTE]: L >= l     :: reason: Appendix-G - 'The GPS interpolation'
l_max = 5       # Number of partial waves: l_max+1            [NOTE]: l_max >= l :: reason: Section-2.3.4 - 'The Partial wave expansion'

x = np.linspace(-1, 1, 199)
r = f(x)

roots, weights = np.polynomial.legendre.leggauss(L+1)       # Gauss-Legendre collocation points (or, nodes) and quadrature weights.
theta = np.linspace(0, 2*np.pi, 200)
R, Theta = np.meshgrid(r, theta)
wavefunction = Psi(l, R, Theta)


gl_array = []
for l_ind in range(0, l_max+1):
    g_l = np.zeros(len(r))
    for i in range(len(g_l)):
        g_k = 0
        for k in range(len(weights)):
            g_k += weights[k] * legendre(l_ind)(roots[k]) * Psi(l, r[i], np.arccos(roots[k]))
        g_l[i] = g_k
    gl_array.append(g_l * np.sqrt(np.pi * (2*l_ind + 1)))

psi = np.zeros((len(theta), len(r)))
for j in range(len(theta)):
    for l_ind in range(len(gl_array)):
        psi[j] += gl_array[l_ind] * Y_l0(l_ind, theta[j])

for l_ind in range(len(gl_array)):
    ax1.plot(r, gl_array[l_ind], 'o-', label=rf'g$_{{{l_ind}}}$(r)')
ax1.plot(r, Psi(l, r, 0) / (np.sqrt((2 * l + 1) / (4 * np.pi))), color='deeppink', label='radial part: ψ(r, θ=0) ', lw=2.5)         # plotting only the radial part at theta=0
ax2.contourf(R * np.sin(Theta), R * np.cos(Theta), wavefunction, 100, cmap='jet')
ax3.contourf(R * np.sin(Theta), R * np.cos(Theta), psi, 100, cmap='jet')

ax2.set_title(r"ψ(r, θ) = $\sqrt{\frac{2\ell + 1}{4\pi}}$ 2r e$^{-r^2}$ P$_\ell$(cosθ); $\ell$=" +f"{l}", pad=30, fontsize=20)
ax3.set_title(fr'ψ(r$_i$, θ$_j$) = $\sum_{{\ell=0}}^{{l_{{max}}}}$g$_{{\ell}}$(r$_i$)P$_{{\ell}}$(cosθ$_j$)', pad=30, fontsize=20)
ax1.legend(loc='upper right', ncol=2, fontsize=15, framealpha=0.5, edgecolor='k')
fig2.subplots_adjust(top=0.86, bottom=0.063)
plt.show()