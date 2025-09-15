"""
File: Partial_wave_expansion_gen_m.py
Project: HHG-SaDAS
Code Description:
    | *** [VERY IMPORTANT CODE] ***
    | This script demonstrates the 'generalised partial wave expansion' of an arbitrary wavefunction.
    |
    | The idea of partial wave expansion is the same as that explained in Partial_wave_Expansion_m=0.py,
    | whereas here, the partial wave expansion is generalised over the m quantum number.
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
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""


import numpy as np
from scipy.special import lpmv
import matplotlib.pyplot as plt
from scipy.special import factorial
import Assistant.Decorate_axes.decorate_axes_D as da

fig1 = plt.figure()
fig2 = plt.figure()
ax1 = fig1.add_subplot(111)             # g_l(r)
ax2 = fig2.add_subplot(121)             # psi(r, θ)
ax3 = fig2.add_subplot(122)             # psi(r, θ) expanded in Legendre Polynomial
da.decorate_2d([ax1])
da.decorate_polar([ax2, ax3])



def f(x):
    L_map = 1
    r_max = 3
    alpha = 2 * L_map / r_max
    return L_map*(1 + x) / (1 - x + alpha)


# NOTE:: I have to include proper starting value of l.
#        Because when -l <= m <= l. means when m=1, l cannot be 0.
#        Example: if m=0, l=(0, L)
#                 if m=1, l=(1, L)

def C_fact(l, m):
    """
    Orthogonality constant factor of P_lm
    """
    return 2 * factorial(l+m) / ((2*l+1)*factorial(l-m))

def N_fact(l, m):
    """
    Normalization constant of Y_lm
    """
    return (-1)**m * np.sqrt((2*l+1) * factorial(l-m) / (4*np.pi * factorial(l+m)))

def a_legendre(l, m, x):
    """
    Associated Legendre polynomial without Condon-Shortley phase
    """
    return lpmv(m, l, x) * (-1)**m

def Y_lm(l, m, x):
    """
    Y_lm(x) = N_lm * P_lm(x)
    [NOTE]: Here phi = 0.
    """
    return N_fact(l, m) * a_legendre(l, m, x)


l = 4
m = 3

def Psi(r, theta):
    return np.exp(-r**2) * Y_lm(l, m, np.cos(theta))


x = np.linspace(-1, 1, 199)
r = f(x)

theta = np.linspace(0, 2*np.pi, 200)
R, Theta = np.meshgrid(r, theta)
wavefunction = Psi(R, Theta)
ax2.contourf(R * np.sin(Theta), R * np.cos(Theta), wavefunction, 100, cmap='jet')

L  = 2           # Number of angular collocation points: L+1   [NOTE]: L >= l       :: reason: Appendix-G - 'The GPS interpolation'
l_max = 2        # Number of partial waves: l_max+1            [NOTE]: l_max >= l-m :: reason: Section-2.3.4 - 'The Partial wave expansion'
roots, weights = np.polynomial.legendre.leggauss(L+1)

gl_array = []
for l_ind in range(0, l_max+1):
    l_eff = l_ind + m
    print(l_eff)
    g_l = np.zeros(len(r))
    for i in range(len(g_l)):
        g_k = 0
        for k in range(len(weights)):
            g_k += weights[k] * a_legendre(l_eff, m, roots[k]) * Psi(r[i], np.arccos(roots[k]))
        g_l[i] = g_k
    gl_array.append(g_l / (N_fact(l_eff, m) * C_fact(l_eff, m)))

psi = np.zeros((len(theta), len(r)))
for j in range(len(theta)):
    for l_ind in range(len(gl_array)):
        l_eff = l_ind + m
        psi[j] += gl_array[l_ind] * Y_lm(l_eff, m, np.cos(theta[j]))

for l_ind in range(len(gl_array)):
    ax1.plot(r, gl_array[l_ind], 'o-', label=rf'g$_{{{l_ind+m}}}$(r)')
ax3.contourf(R * np.sin(Theta), R * np.cos(Theta), psi, 100, cmap='jet')


ax2.set_title(fr'ψ(r, θ) = e$^{{-r^2}} \cdot P_{{\ell m}}$(cosθ); l, m = {l}, {m}', pad=40, fontsize=20)
ax3.set_title(fr'ψ(r$_i$, θ$_j$) = $\sum_{{\ell=m}}^{{l_{{max}}+m}}$g$_{{\ell}}$(r$_i$)P$_{{\ell m}}$(cosθ$_j$)', pad=40, fontsize=20)
ax1.legend(loc='upper right', ncol=2, fontsize=15, framealpha=0.5, edgecolor='k')
fig2.subplots_adjust(top=0.86, bottom=0.063)
plt.show()
