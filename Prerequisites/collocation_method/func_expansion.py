"""
File: func_expansion.py
Project: HHG-SaDAS

Code Description:
    | This code is aimed at providing a better understanding of:
    |   • Gauss–Lobatto cardinal functions
    |   • Pseudo-spectral interpolation
    |   • Derivative matrices
    |
    | It verifies the correctness of the computed first- and second-order
    | derivative matrices, d(1)_ij and d(2)_ij.


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS

--------------------------------------------------------------------------------
Notes:
- This is the first effort to calculate the collocation points.
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

def find_roots(f, x):
    """
    :param f: function; which returns the value of function.
    :param x: array; passed as argument of f.
    :return: I) array; All 0 roots of 'f' function, between min(x) & max(x)
             II) array; The values of 'f' function at those roots.
    """
    pks_at = find_peaks(-f(x) ** 2)[0]
    rts = fsolve(f, x[pks_at])
    val_at_rt = f(rts)
    return rts, val_at_rt

fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
fig4 = plt.figure()
ax1 = fig1.add_subplot(111)              # Roots of P'(N, x)
ax2 = fig2.add_subplot(211)              # Collocation method and interpolation
ax3 = fig2.add_subplot(212)              # Error between exact and interpolated function
ax4 = fig3.add_subplot(111)              # cardinal function: gj(x)
ax5 = fig4.add_subplot(211)              # First derivative of the cardinal function: gj'(x)
ax6 = fig4.add_subplot(212)              # Second derivatives of the cardinal function: gj"(x)
da.decorate_2d([ax1, ax2, ax3, ax4, ax5, ax6])



x, dx = np.linspace(-1, 1, 250, retstep=True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~: Getting P(l, x), P'(l, x) & the collocation points:~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = 20
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

def a_derivative_P_N_single(x):
    return N*(N+1) * (legendre(N-1)(x) - legendre(N+1)(x)) / (1-x**2) / (2*N+1)

def root_func(x):
    return legendre(N-1)(x) - legendre(N+1)(x)

PN_poly = P_N(x)
rts_PN_poly = find_roots(P_N, x)[0]


root_x = np.linspace(-1, 1, 2000)
root_PN_deriv_array = root_func(root_x)
a_PN_deriv_array = a_derivative_P_N(x)
pks_at = find_peaks(-root_PN_deriv_array ** 2)[0]
colloc_pt = np.array([fsolve(root_func, root_x[pks_at[j]])[0] for j in range(len(pks_at))])
# colloc_pt = np.insert(colloc_pt, [0, -1], [-1, 1])

colloc_nop = len(colloc_pt)
print('Order of polynomial (N)      :', N)
print('number of collocation points :', colloc_nop)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~: Calculating gj(x) and collocation expansion :~~~~~~~~~~~~~~~~~~~~~~~~~~~
def gj(xj):
    """
    Cardinal function
    xj :  Collocation points [the roots of P'(N,x)].
    """
    # if xj == -1:
    #     ret1 = np.zeros(len(x)); ret1[0]=1
    #     return ret1
    # if xj == 1:
    #     ret2 = np.zeros(len(x)); ret2[-1]=1
    #     return ret2
    # else:
    term1 = 1 / (N * (N+1))
    term2 = (1-x**2) / (xj - x)
    term3 = a_PN_deriv_array / P_N(xj)
    return term1 * term2 * term3

def phi(x):
    return np.sin(2*np.pi*x)
    # return np.sin(5*x)

phi_N = np.zeros(len(x))
cardinal_func_arr = []
for j in range(colloc_nop):
    cardinal_func = gj(colloc_pt[j])
    phi_N += cardinal_func * phi(colloc_pt[j])
    cardinal_func_arr.append(cardinal_func)
    ax4.plot(x, cardinal_func, label=rf'g$_{{{j}}}$(x)')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~: Checking the first derivative :~~~~~~~~~~~~~~~~~~~~~~~~~~~
def d1(i, j):
    # paper in support     : https://journals.aps.org/pra/pdf/10.1103/PhysRevA.59.2864    (D.A.Telnov & Shih-I Chu, Phys.Rev.A59,2864(1999))
    #                      : https://journals.aps.org/pra/pdf/10.1103/PhysRevA.76.043412  (D.A.Telnov & Shih-I Chu, Phys.Rev.A76,043412(2007))
    # paper not in support : https://journals.aps.org/pra/pdf/10.1103/PhysRevA.91.063412  (Bohmian trajectory)
    #                      : https://doi.org/10.1002/qua.26245                            (Review GPSM)
    if i != j:
        return 1 / (colloc_pt[i] - colloc_pt[j])
    # if i == j == 0:
    #     return -(N + 1)*N / 4
    # if i == j == N:
    #     return (N + 1)*N / 4
    else:
        return 0                    # if none of the above conditions are true


def gj_p(i, j):
    """
    First derivative of cardinal function
    """
    return d1(i, j) * P_N(colloc_pt[i]) / P_N(colloc_pt[j])


# which cardinal function's derivative do you want to see...
j = 5               # has to be: 0 < j < N-1
gj_p_exact = np.gradient(cardinal_func_arr[j], x)   # less sophisticated derivative
gj_p_matrix = [gj_p(i, j) for i in range(colloc_nop)]


xt = colloc_pt[j] - 0.0001          # added 0.0001 (+ or -) so that term1 won't just blow up. Majorly to act as limit: x -> xj
term1 = -2 * P_N(xt) / ((xt - colloc_pt[j])**2 * P_N(colloc_pt[j]))
term2 = (1 / (N * (N + 1))) * ((2 * (xt**2 - 1) / ((xt - colloc_pt[j])**3)) +
                               (N * (N + 1) / (xt - colloc_pt[j]))) * (a_derivative_P_N_single(xt) / P_N(colloc_pt[j]))
gj_dp_test = term1 + term2                                # double derivative at x=xj
gj_dp_orig = -N*(N+1) / (3*(1 - colloc_pt[j]**2))

print("original gj''(xj)      :", gj_dp_orig)
print("test gj''(xj)          :", gj_dp_test)
print('REMARK: this implies the analytically calculated second derivative is correct. SORTED: DONE')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~: Checking the second derivative :~~~~~~~~~~~~~~~~~~~~~~~~~~~
def d2(i, j):
    if i != j:
        return -2 / (colloc_pt[i] - colloc_pt[j])**2
    else:
        # return N*(N+1) / (1 - colloc_pt[i]**2)
        return -N*(N+1) / (3*(1 - colloc_pt[i]**2))

def gj_pp(i, j):
    """
    Second derivative of cardinal function
    """
    return d2(i, j) * P_N(colloc_pt[i]) / P_N(colloc_pt[j])

gj_pp_exact = np.gradient(gj_p_exact, x)
gj_pp_matrix = [gj_pp(i, j) for i in range(colloc_nop)]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~: Plotting :~~~~~~~~~~~~~~~~~~~~~~~~~~~
ax1.plot(x, a_PN_deriv_array, label=f"P'(N={N}, x)")
ax1.scatter(colloc_pt, np.zeros(len(colloc_pt)), color='red', zorder=2)
ax2.plot(x, phi(x), 'o-', markersize=4, color='deeppink', label='φ(x)')
ax2.plot(x, phi_N, 'o-', markersize=2, label=r'φ$_N$(x)')
ax3.plot(x, phi(x) - phi_N, label=r'φ(x) - φ$_N$(x)')
ax5.plot(x, cardinal_func_arr[j] * max(gj_p_exact) / max(cardinal_func_arr[j]), label=rf"g$_{{{j}}}$(x)")
ax5.plot(x, gj_p_exact, label=rf"g$_{{{j}}}$'(x)")
ax5.scatter(colloc_pt, gj_p_matrix, s=60, color='yellow', zorder=2)
ax6.plot(x, cardinal_func_arr[j] * max(gj_pp_exact) / max(cardinal_func_arr[j]), label=rf"g$_{{{j}}}$(x)")
ax6.plot(x, gj_pp_exact, label=rf"g$_{{{j}}}$''(x)")
ax6.scatter(colloc_pt, gj_pp_matrix, s=50, color='yellow', zorder=2)

ax1.set_ylim(-20, 20)
ax4.set_ylim(-0.4, 1.6)
ax4.scatter(colloc_pt, np.zeros(len(colloc_pt)), s=45, color='red', zorder=2)
ax3.scatter(colloc_pt, np.zeros(len(colloc_pt)), s=45, color='red', zorder=2)
ax5.scatter(colloc_pt, np.zeros(len(colloc_pt)), s=45, color='red', zorder=2)
ax6.scatter(colloc_pt, np.zeros(len(colloc_pt)), s=45, color='red', zorder=2)

ax1.legend(loc='upper center', fontsize=12, framealpha=0.5, edgecolor='k')
ax2.legend(loc='upper right', fontsize=12, framealpha=0.5, edgecolor='k')
ax3.legend(loc='upper right', fontsize=12, framealpha=0.5, edgecolor='k')
ax4.legend(loc='upper right', ncol=4, fontsize=12, framealpha=0.5, edgecolor='k')
ax5.legend(loc='upper right', fontsize=12, framealpha=0.5, edgecolor='k')
ax6.legend(loc='upper right', fontsize=12, framealpha=0.5, edgecolor='k')
plt.show()