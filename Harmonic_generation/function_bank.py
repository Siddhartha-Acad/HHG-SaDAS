"""
File: function_bank.py
Project: HHG-SaDAS
Code Description:
    | This contains all functions that are used in HHG-SaDAS package.
    | All other simulating codes, fetch functions from this script.


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

from scipy.special import legendre
from Harmonic_generation.parameters import *


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               radial mapping function f(x)               |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def f(x, Lmap=L_map):
    """
    Nonlinear radial mapping function.

    :param x: Input variable, typically defined in the range [-1, 1] (float or array-like).
    :param Lmap: Mapping parameter that controls the scaling of the transformation (float).
    :return: Transformed value(s) according to the nonlinear radial mapping (float or numpy.ndarray).

    References
    ----------
    For details, see Section 2.2.2 -- 'Grid discretization and nonlinear mapping'
    """
    alpha = 2 * Lmap / r_max
    return Lmap * (1 + x) / (1 - x + alpha)


def f_p(x, Lmap=L_map):
    """
    First derivative of the nonlinear radial mapping function.

    :param x: Input variable, typically defined in the range [-1, 1] (float or array-like).
    :param Lmap: Mapping parameter that controls the scaling of the transformation (float).
    :return: Value(s) of the derivative of the nonlinear radial mapping (float or numpy.ndarray).
    """
    alpha = 2 * Lmap / r_max
    return Lmap * (alpha + 2) / (1 - x + alpha)**2


def P_N(x):
    """
    Legendre polynomial of degree N at x.

    :param x: Point(s) of evaluation (float or array-like).
    :param N: Polynomial degree (int).
    :return: Value(s) of P_N(x) (float or numpy.ndarray).
    """
    return legendre(N)(x)


def d2(i, j):
    """
    Second-order derivative matrix of Gauss–Lobatto cardinal functions.
    This function returns the (i, j) element of the second derivative matrix
    associated with Gauss–Lobatto collocation points.

    :param i: Row index (int).
    :param j: Column index (int).
    :return: Value of the (i, j) entry of the second derivative matrix (float).

    References
    ----------
    For details, see Section 2.2.7 -- 'Applying GPSM to construct the matrix Hamiltonian'
    """
    if i != j:
        return -2 / (colloc_pt[i] - colloc_pt[j])**2
    else:
        return -N*(N+1) / (3*(1 - colloc_pt[i]**2))


def V_eff(l, x):
    """
    Effective potential for a particle in a central field.

    :param l: Angular momentum quantum number (int).
    :param x: Radial coordinate (float).
    :return: Effective potential value at x (float).
    """
    return -1 / x + l*(l+1) / (2*x**2)


def H(l, i, j, model='SAE-M2'):
    term1 = -0.5 * (1 / f_p(colloc_pt[i])) * d2(i, j) * (1 / f_p(colloc_pt[j]))
    if i != j:
        return term1
    if model == 'SAE-M1':
        term2 = l * (l + 1) / (2 * f(colloc_pt[i]) ** 2) + potential_V_SAE_M1(f(colloc_pt[i]), atom=evolving_atom)
        return term1 + term2
    elif model == 'SAE-M2':
        term2 = l * (l + 1) / (2 * f(colloc_pt[i]) ** 2) + potential_V_SAE_M2(f(colloc_pt[i]), atom=evolving_atom)
        return term1 + term2

def S(E_l, A_l, i, j):
    return sum(A_l[k][i] * A_l[k][j] * np.exp(-1j * E_l[k] * dt / 2) for k in range(len(E_l)))

def Up(E0_au, w0):                  # Ponderomotive energy
    return E0_au**2 / (4 * w0**2)

def Keldysh(Ip_au, Up_au):          # Keldysh Parameter
    return np.sqrt(Ip_au / (2*Up_au))

def N_cutoff(Ip_au, Up_au):         # Cut-off Harmonic
    return (Ip_au + 3.17*Up_au) / w0

def generate_states(l):
    orbital_types = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i', 7: 'j', 8: 'k', 9: 'l', 10: 'm'}
    return [str(i) + orbital_types.get(l, '') for i in range(l + 1, 200)]


def E_field(t):
    return E0_au * np.sin(w0*t) * (np.sin(w0*t / (2*cpp))) ** 2

def V_int(r, theta, t):
    return -E_field(t) * r * np.cos(theta)

def Absorber_func(r):
    if 0 < r <= r0: return 1
    elif r0 < r < r_max:
        return np.cos(np.pi * (r - r0) / (2 * (r_max - r0))) ** 0.25


def potential_V_SAE_M1(r, atom='Ne'):
    params = atomic_params_SAE_M1.get(atom)
    if params is None: raise ValueError(f"Atom or ion '{atom}' not found in table.")
    Zc, a1, a2, a3, a4, a5, a6 = params["Zc"], params["a1"], params["a2"], params["a3"], params["a4"], params["a5"], params["a6"]
    return -(Zc + a1*np.exp(-a2*r) + a3*r*np.exp(-a4*r) + a5*np.exp(-a6*r)) / r

def potential_V_SAE_M2(r, atom='Ne'):
    params = atomic_params_SAE_M2.get(atom)
    if params is None: raise ValueError(f"Atom or ion '{atom}' not found in table.")

    C0, Zc, c = params['C0'], params['Zc'], params['c']
    a1, a2, a3 = params['a1'], params['a2'], params['a3']
    b1, b2, b3 = params['b1'], params['b2'], params['b3']

    V_long = -C0 / r
    V_short = -Zc * np.exp(-c * r) / r
    V_shell = (a1 * np.exp(-b1 * r)) + (a2 * np.exp(-b2 * r)) + (a3 * np.exp(-b3 * r))
    return V_long + V_short - V_shell       # overall -ve sign in V_shell.



def dydx(integrand, x):
    """
    FIRST ORDER DERIVATIVE

    :param integrand: Values of the function at discrete points (array-like).
    :param x: Corresponding x-values (array-like).
    :return: Derivative of the function (numpy.ndarray).
    """
    nop_x = len(x)
    dy_dx = np.zeros_like(integrand, dtype=np.float64)
    dy_dx[0] = (integrand[1] - integrand[0]) / (x[1] - x[0])
    for i in range(1, nop_x - 1):
        dy_dx[i] = (integrand[i + 1] - integrand[i - 1]) / (x[i + 1] - x[i - 1])
    dy_dx[-1] = (integrand[-1] - integrand[-2]) / (x[-1] - x[-2])
    return dy_dx