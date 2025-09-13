"""
File: GPSM_interp_matrix.py
Project: HHG-SaDAS

Code Description:
    | This code is aimed at providing a better understanding of:
    |   • Gauss–Lobatto cardinal functions
    |   • Pseudo-spectral interpolation
    |   • Derivative matrices
    |   • Interpolation matrix
    |
    | It verifies the correctness of the computed first- and second-order
    | derivative matrices, d(1)_ij and d(2)_ij.


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS

--------------------------------------------------------------------------------
Notes:
- It calculates and shows the GPSM interpolation matrix.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""


import numpy as np
from Assistant.Decorate_axes import decorate_axes_D as da
import matplotlib.pyplot as plt
from scipy.special import legendre
from scipy.signal import find_peaks
from scipy.optimize import fsolve

fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
ax1 = fig1.add_subplot(111)              # Interpolation matrix (notebook date: 08/02/2025)
ax2 = fig2.add_subplot(211)              # Collocation method and interpolation
ax3 = fig2.add_subplot(212)              # Error between exact and interpolated function
ax4 = fig3.add_subplot(111)              # cardinal function: gj(x)
da.decorate_2d([ax1, ax2, ax3, ax4])

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


x = np.linspace(-1, 1, 50)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~: Getting P(l, x), P'(l, x) & the collocation points:~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = 50
root_x = np.linspace(-1, 1, 2000)
root_PN_deriv_array = root_func(root_x)
a_PN_deriv_array = a_derivative_P_N(x)
pks_at = find_peaks(-root_PN_deriv_array ** 2)[0]
colloc_pt = np.array([fsolve(root_func, root_x[pks_at[j]])[0] for j in range(len(pks_at))])

colloc_nop = len(colloc_pt)
print('Order of polynomial (N)      :', N)
print('number of collocation points :', colloc_nop)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~: Calculating gj(x) and collocation expansion :~~~~~~~~~~~~~~~~~~~~~~~~~~~
def gj(xj):
    """
    Cardinal function
    xj :  Collocation points [the roots of P'(N,x)].
    """
    term1 = 1 / (N * (N+1))
    term2 = (1-x**2) / (xj - x)
    term3 = a_PN_deriv_array / P_N(xj)
    return term1 * term2 * term3

def phi(x):
    return np.sin(2*np.pi*x)
    # return np.sin(5*x)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~: Plotting :~~~~~~~~~~~~~~~~~~~~~~~~~~~
phi_N = np.zeros(len(x))
cardinal_func_matrix = []
for j in range(colloc_nop):
    cardinal_func = gj(colloc_pt[j])
    phi_N += cardinal_func * phi(colloc_pt[j])
    cardinal_func_matrix.append(cardinal_func)
    ax4.plot(x, cardinal_func, label=rf'g$_{{{j}}}$(x)')
cardinal_func_matrix = np.array(cardinal_func_matrix).T     # each vertical column represents each cardianl function

ax1.pcolormesh(cardinal_func_matrix)
ax1.set_aspect('equal')

ax2.plot(x, phi(x), 'o-', markersize=4, color='deeppink', label='φ(x)')
ax2.plot(x, phi_N, 'o-', markersize=2, label=r'φ$_N$(x)')
ax3.plot(x, phi(x) - phi_N, label=r'φ(x) - φ$_N$(x)')

ax4.set_ylim(-0.4, 1.6)
ax4.scatter(colloc_pt, np.zeros(len(colloc_pt)), s=45, color='red', zorder=2)
ax3.scatter(colloc_pt, np.zeros(len(colloc_pt)), s=45, color='red', zorder=2)

ax2.legend(loc='upper right', fontsize=12, framealpha=0.5, edgecolor='k')
ax3.legend(loc='upper right', fontsize=12, framealpha=0.5, edgecolor='k')
plt.show()