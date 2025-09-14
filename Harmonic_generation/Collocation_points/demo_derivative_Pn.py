"""
File: demo_derivative_Pn.py
Project: HHG-SaDAS
Code Description:
    | Demonstrates the basic principles of Algorithm-3 (Algo-3) for analyzing collocation points.
    | The script:
    |   - Maps a dense grid of x-values using a reversed mapping function.
    |   - Computes the numerical derivative of P_N(x).
    |   - Computes the analytical derivative using the closed-form expression.
    |   - Plots comparisons between numerical and analytical derivatives,
    |     the mapping function, and the root function P_{N-1}(x) - P_{N+1}(x).

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- The script is for demonstration purposes only and illustrates the base
  methodology of Algo-3, including error visualization and mapping strategies.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.special import legendre
from Assistant.Decorate_axes import decorate_axes_D as da

fig1 = plt.figure()
ax1 = fig1.add_subplot(231)
ax2 = fig1.add_subplot(234)
ax3 = fig1.add_subplot(232)
ax4 = fig1.add_subplot(235)
ax5 = fig1.add_subplot(233)
ax6 = fig1.add_subplot(236)
da.decorate_2d([ax1, ax2, ax3, ax4, ax5, ax6])

def f_rev(x_array):             # f(x) reversed to have dense grid towards far.
    r_max = 1; L_map = 0.2
    alpha = 2 * L_map / r_max
    map_func = L_map*(1 + x_array) / (1 - x_array + alpha)
    return -map_func[::-1] + r_max

def P_N(x):
    return legendre(N)(x)

def a_derivative_P_N(x):
    if x[0] == -1 and x[-1] == 1:
        x = x[1: -1]
        ret = N*(N+1) * (legendre(N-1)(x) - legendre(N+1)(x)) / (1-x**2) / (2*N+1)
        return np.concatenate(([N*(N+1)/2 * (-1)**(N-1)], ret, [N*(N+1)/2]))
    if x[0] == 0 and x[-1] == 1:
        x = x[0: -1]
        ret =  N*(N+1) * (legendre(N-1)(x) - legendre(N+1)(x)) / (1-x**2) / (2*N+1)
        return np.append(ret, N * (N + 1) / 2)

def root_func(x):
    return legendre(N-1)(x) - legendre(N+1)(x)

N = 20
x = np.linspace(-1, 1, 500); dx = x[1] - x[0]
x_mapped = f_rev(x)

n_deriv = np.array([derivative(P_N, float(x[i]), float(dx), n=1) for i in range(len(x))])
a_deriv = a_derivative_P_N(x)


ax1.plot(x, n_deriv, 'o-', markersize=10, label=r"Numerical P$^{'}_{N}$(x)")
ax1.plot(x, a_deriv, 'o-', color='m', label=r"Analytical P${'}_{N}$(x)")
ax1.plot(x, 20*root_func(x), 'o-', markersize=3, color='orange', label=r'P$_{N-1}$(x) - P$_{N+1}$(x)')
ax2.plot(x, x_mapped, 'o-', label='mapped x')
ax3.plot(x[len(x) // 2:], a_deriv[len(x) // 2:], 'o-', color='m', label=r"Analytical P${'}_{N}$(x)")
ax4.plot(x_mapped, a_derivative_P_N(x_mapped), 'o-', color='#83C167', label=r"Analytical P${'}_{N}$(f(x))")
ax5.plot(x[len(x) // 2:], root_func(x)[len(x) // 2:], 'o-', color='orange', label=r'P$_{N-1}$(x) - P$_{N+1}$(x)')
ax6.plot(x_mapped, root_func(x_mapped), 'o-', color='crimson', label=r'P$_{N-1}$(f(x)) - P$_{N+1}$(f(x))')

da.da_legend([ax1, ax2, ax3, ax4, ax5, ax6], loc='upper center', fontsize=11)

ax1_ylim = ax1.get_ylim(); ax2_ylim = ax2.get_ylim()
ax3_ylim = ax3.get_ylim(); ax4_ylim = ax4.get_ylim()
ax5_ylim = ax5.get_ylim(); ax6_ylim = ax6.get_ylim()

# ax1.set_ylim(ax1_ylim[0] * 1.5, ax1_ylim[1] * 1.5)
ax1.set_ylim(-17, 17)
ax5.set_ylim(ax5_ylim[0] * 1.5, ax5_ylim[1] * 1.5)
ax6.set_ylim(ax6_ylim[0] * 1.5, ax6_ylim[1] * 1.5)
fig1.subplots_adjust(top=0.954, bottom=0.05, left=0.034, right=0.983, wspace=0.133, hspace=0.145)

plt.show()