"""
File: parameters_and_functions.py
Project: HHG-SaDAS
Code Description:
    | This contains all parameters that defines the entire system and numerical requirements.
    | All other simulating codes, fetch parameters and functions from this script.


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

import warnings
import numpy as np
from pathlib import Path
from scipy.special import legendre
from Atomic_units import Int_0, omega_au, a0, T0
this_dir = Path(__file__).resolve().parent              # Relative path system

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
print('*** RuntimeWarning     : Blocked from parameters_and_functions.py ***')
print('*** DeprecationWarning : Blocked from parameters_and_functions.py ***\n')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Atom, SAE and Confinement            |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n = 2; l = 0; m = 0         # defines initial state. [NOTE]: for a given l, n always starts from 1. Ex: 1s=(1, 0, 0); 2pz=(1, 1, 0); 4px=(3, 1, 1)
evolving_atom = 'He'        # Atoms are listed down in 'SAE dataset' section.
SAE_model = 'SAE-M1'        # Single active electron model; option: SAE_model = 'SAE-M1' or 'SAE-M2'. [NOTE]: For 'Xe' always use 'SAE-M1'

confined = False                    # whether the atom is confined or not?
confinement_model = 'P-Gau'         # which type of confinement potential?
conf_info_string = 'till_empty'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                GPSM Parameters                 |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = 200                     # P'_N(xj) = 0 ; radial grid size: len(colloc_pt) = N-1
L = 20                      # must be >= l ; number of S matrix gl(r) in partial wave expansion.
l_max = 20                  # Number of partial waves = number of S-matrices = l_max+1
k_max = 50                  # number of GPSM states (maximum k index) in S matrix
L_map = 20; r_max = 200     # radial mapping parameters
r0 = 150                    # absorber layer thickness: (r_max - r0) a.u.

colloc_file = f'Algo-3_N={N}_AnaDeriv_collocation_points.txt'
colloc_file = this_dir.parent / 'Harmonic_Generation' / 'Collocation_points' / 'AnaDeriv_Colloc_pt' / colloc_file
colloc_pt = np.loadtxt(colloc_file, skiprows=1, usecols=0)

int_w = 2 / (N * (N + 1) * (legendre(N)(colloc_pt))**2)           # Gauss-Lobatto  Quadrature weights: w_j
roots, weights = np.polynomial.legendre.leggauss(L+1)             # Gauss-Legendre Quadrature weights and collocation points (or, nodes): x_k
theta_k = np.arccos(roots)                                        # Angular collocation points: cos(theta_k) = x_k


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          LASER and temporal grid info          |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lambda_nm = 1064                                # wavelength (nm)
I0 = 5 * 10**13                                 # Intensity (W/cm2)
I0_au = I0 / Int_0                              # Intensity (a.u)
E0_au = np.sqrt(I0_au)                          # Field intensity (a.u)
w0 = omega_au(lambda_nm); T = 2 * np.pi / w0    # Angular frequency and time period.
cpp = 60; tf = cpp*T; dt = 0.1                  # cpp = cycles per pulse.
t = np.arange(0, tf+dt, dt)                     # total number of time steps.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  SAE dataset                   |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atomic_params_SAE_M1 = {        # Ref: X. M. Tong and C. D. Lin, J. Phys. B: At. Mol. Opt. Phys., 38, 2593 (2005).
    "H"  :  {"Zc": 1.0, "a1": 0.000,  "a2": 0.000,  "a3": 0.000,   "a4": 0.000, "a5": 0.000,  "a6": 0.000},
    "He" :  {"Zc": 1.0, "a1": 1.231,  "a2": 0.662,  "a3": -1.325,  "a4": 1.236, "a5": -0.231, "a6": 0.480},
    "Ne" :  {"Zc": 1.0, "a1": 8.069,  "a2": 2.148,  "a3": -3.570,  "a4": 1.986, "a5": 0.931,  "a6": 0.602},
    "Ar" :  {"Zc": 1.0, "a1": 16.039, "a2": 2.007,  "a3": -25.543, "a4": 4.525, "a5": 0.961,  "a6": 0.443},
    "Xe" :  {"Zc": 1.0, "a1": 51.356, "a2": 2.112,  "a3": -99.927, "a4": 3.737, "a5": 1.644,  "a6": 0.431},
    # "Rb" :  {"Zc": 1.0, "a1": 24.023, "a2": 11.107, "a3": 115.200, "a4": 6.629, "a5": 11.977, "a6": 1.245},
    # "Ne+":  {"Zc": 2.0, "a1": 8.043,  "a2": 2.715,  "a3": 0.506,   "a4": 0.982, "a5": -0.043, "a6": 0.401},
    # "Ar+":  {"Zc": 2.0, "a1": 14.989, "a2": 2.217,  "a3": -23.606, "a4": 4.585, "a5": 1.011,  "a6": 0.551}
}

atomic_params_SAE_M2 = {        # Ref: R. Reiff, T. Joyce, A. Jaroń-Becker, and A. Becker, J. Phys. Commun., 4, 065011 (2020).
    "H"   : {"C0": 1, "Zc": 0,  "c": 0.000,   "a1": 0.000,    "a2": 0.000,    "a3": 0.000,  "b1": 0.000,   "b2": 0.000,   "b3": 0.000},
    "He"  : {"C0": 1, "Zc": 1,  "c": 2.0329,  "a1": 0.3953,   "a2": 0.000,    "a3": 0.000,  "b1": 6.1805,  "b2": 0.000,   "b3": 0.000},
    "Ne"  : {"C0": 1, "Zc": 9,  "c": 0.8870,  "a1": -9.9286,  "a2": -5.9950,  "a3": 0.000,  "b1": 1.3746,  "b2": 3.7963,  "b3": 0.000},
    "Ar"  : {"C0": 1, "Zc": 17, "c": 0.8103,  "a1": -15.9583, "a2": -27.7467, "a3": 2.1768, "b1": 1.2305,  "b2": 4.3946,  "b3": 86.7179},
    # "Li"  : {"C0": 1, "Zc": 2,  "c": 15.9594, "a1": 9.1124,   "a2": 19.3145,  "a3": 0.000,  "b1": 3.6040,  "b2": 11.3082, "b3": 0.000},
    # "Be"  : {"C0": 1, "Zc": 3,  "c": 2.0481,  "a1": 0.5294,   "a2": 0.3219,   "a3": 0.000,  "b1": 0.8475,  "b2": 37.5567, "b3": 0.000},
    # "Na"  : {"C0": 1, "Zc": 10, "c": 1.4927,  "a1": -11.3552, "a2": -2.0302,  "a3": 1.6028, "b1": 2.5597,  "b2": 10.1463, "b3": 47.9555},
    # "Mg"  : {"C0": 1, "Zc": 11, "c": 1.4248,  "a1": -14.5892, "a2": -1.9433,  "a3": 1.8141, "b1": 2.7001,  "b2": 12.3150, "b3": 51.7100},
    # "Ar+" : {"C0": 2, "Zc": 16, "c": 0.8698,  "a1": -16.0391, "a2": -26.9860, "a3": 2.1780, "b1": 1.3146,  "b2": 4.4514,  "b3": 88.3315},
    # "Ar2+": {"C0": 3, "Zc": 15, "c": 0.8792,  "a1": -16.4007, "a2": -26.6805, "a3": 2.1681, "b1": 1.3486,  "b2": 4.4656,  "b3": 90.3068},
    # "Ar3+": {"C0": 4, "Zc": 14, "c": 0.9445,  "a1": -16.4800, "a2": -25.8243, "a3": 2.1550, "b1": 1.4521,  "b2": 4.5171,  "b3": 94.5151},
    # "Ar4+": {"C0": 5, "Zc": 13, "c": 0.8529,  "a1": -17.4441, "a2": -26.0893, "a3": 2.1135, "b1": 1.4141,  "b2": 4.4613,  "b3": 101.5018},
    # "Ar5+": {"C0": 6, "Zc": 12, "c": 0.8929,  "a1": -17.5407, "a2": -25.4398, "a3": 2.0818, "b1": 1.5024,  "b2": 4.4823,  "b3": 108.4695}
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Function Bank :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           radial mapping function f(x) & f'(x)           |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def f(x, Lmap=L_map):
    r"""
    Nonlinear radial mapping function.

    .. math::
        f(x) = \frac{L_{\mathrm{map}} \, (1 + x)}{\,1 - x + \alpha\,},
        \qquad \alpha = \frac{2 L_{\mathrm{map}}}{r_{\max}}

    :param x: Input variable, typically in the range [-1, 1] (float or array-like).
    :param Lmap: Mapping parameter controlling the scaling of the transformation (float).
    :return: Transformed value(s) according to the nonlinear radial mapping (float or ndarray).

    Reference
    ----------
    - Section 2.2.2 — *Grid discretization and nonlinear mapping*
    """
    alpha = 2 * Lmap / r_max
    return Lmap * (1 + x) / (1 - x + alpha)


def f_p(x, Lmap=L_map):
    r"""
    First derivative of the nonlinear radial mapping function.

    .. math::
        f'(x) = \frac{L_{\mathrm{map}} \, (\alpha + 2)}{(1 - x + \alpha)^2},
        \qquad \alpha = \frac{2 L_{\mathrm{map}}}{r_{\max}}

    :param x: Input variable, typically in the range [-1, 1] (float or array-like).
    :param Lmap: Mapping parameter controlling the scaling of the transformation (float).
    :return: Value(s) of the derivative of the nonlinear radial mapping (float or ndarray).
    """
    alpha = 2 * Lmap / r_max
    return Lmap * (alpha + 2) / (1 - x + alpha)**2



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          Second-order derivative matrix : d2_ij          |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def d2(i, j):
    """
    Second-order derivative matrix of Gauss–Lobatto cardinal functions.
    This function returns the (i, j) element of the second derivative matrix
    associated with Gauss–Lobatto collocation points.

    :param i: Row index (int).
    :param j: Column index (int).
    :return: Value of the (i, j) entry of the second derivative matrix (float).

    Reference
    ----------
    - Section 2.2.7 -- *Applying GPSM to construct the matrix Hamiltonian*
    """
    if i != j:
        return -2 / (colloc_pt[i] - colloc_pt[j])**2
    else:
        return -N*(N+1) / (3*(1 - colloc_pt[i]**2))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   H-matrix & S-matrix                    |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def H(l, i, j, model='SAE-M2'):
    r"""
    Real symmetric Hamiltonian matrix element in the radial mapped discrete Gauss-Lobatto collocation grid.

    .. math::
        \left[H^{\ell}\right]_{ij} =
        -\frac{1}{2} \, \frac{1}{f'(x_i)} \, d^{2}_{ij} \, \frac{1}{f'(x_j)}
        \;+\; \delta_{ij} \left[
            \frac{\ell(\ell+1)}{2 f(x_i)^2}
            + V_{\text{SAE}}(f(x_i))
        \right]

    where:
        - :math:`f(x)` is the nonlinear radial mapping function,
        - :math:`d^{2}_{ij}` is the second derivative matrix element,
        - :math:`V_{\text{SAE}}` is the single-active-electron (SAE) potential,
          with model choice ``SAE-M1`` or ``SAE-M2``.

    :param l: Angular momentum quantum number :math:`\ell` (int).
    :param i: Basis index (row) (int).
    :param j: Basis index (column) (int).
    :param model: SAE model to use, either ``'SAE-M1'`` or ``'SAE-M2'`` (str, default='SAE-M2').
    :return: Hamiltonian matrix element :math:`\left[H^{\ell}\right]_{ij}` (float, a.u.).

    References
    ----------
    - Section 2.2 -- *Numerical solution of hydrogen atom*
    - Section 2.2.8 -- *Symmetrization of the Hamiltonian matrix*
    """
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
    r"""
    The S-matrix elements.

    .. math::
        S_{\alpha\beta}(\ell) =
        \langle x_\alpha \,|\, \exp\!\left(-i \hat{h}^{(0)}_\ell(x) \, \delta t / 2\right) \,|\, x_\beta \rangle

    :param E_l: Eigenenergies (array-like, length k_max).
                 Here k_max = len(E_l).
    :param A_l: Eigenvectors (2D array-like).
    :param i: Basis index (int).
    :param j: Basis index (int).
    :return: S-matrix element S[i, j] (complex).

    References
    ----------
    - Appendix C — *Derivation and Consistency of S-matrix formalism*
    - Section 2.3.5 — *Matrix time evolution operator: The S-matrix*
    """
    return sum(A_l[k][i] * A_l[k][j] * np.exp(-1j * E_l[k] * dt / 2) for k in range(len(E_l)))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ponderomotive force;  Keldysh parameter; Harmonic cut-off |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Up(E0_au, w0):
    r"""
    Ponderomotive energy :math:`U_p`.

    .. math::
        U_p = \frac{E_0^2}{4 \omega_0^2}

    :param E0_au: Electric field amplitude (a.u.).
    :param w0: Angular frequency of the laser field (a.u.).
    :return: Ponderomotive energy (a.u.).
    """
    return E0_au**2 / (4 * w0**2)


def Keldysh(Ip_au, Up_au):
    r"""
    Keldysh parameter :math:`\gamma`.

    .. math::
        \gamma = \sqrt{\frac{I_p}{2 U_p}}

    :param Ip_au: Ionization potential (a.u.).
    :param Up_au: Ponderomotive energy (a.u.).
    :return: Keldysh parameter (dimensionless).
    """
    return np.sqrt(Ip_au / (2 * Up_au))


def N_cutoff(Ip_au, Up_au):
    r"""
    Harmonic cut-off order :math:`N_c`.

    .. math::
        N_c = \frac{I_p + 3.17 U_p}{\omega_0}

    :param Ip_au: Ionization potential (a.u.).
    :param Up_au: Ponderomotive energy (a.u.).
    :param w0: Angular frequency of the laser field (a.u.).
    :return: Cut-off harmonic order (dimensionless).
    """
    return (Ip_au + 3.17 * Up_au) / w0



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      LASER electric field and interaction-potential      |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def E_field(t):
    r"""
    Laser electric field with a sine-squared envelope.

    .. math::
        E(t) = E_0 \, \sin(\omega_0 t) \,
        \left[ \sin\!\left(\frac{\omega_0 t}{2 N_c}\right) \right]^2

    :param t: Time (a.u.).
    :return: Electric field amplitude at time t (a.u.).

    References
    ----------
    - Section 2.3 — *Time evolution of the atomic wavefunction interacting with an external strong-field laser*
    """
    return E0_au * np.sin(w0 * t) * (np.sin(w0 * t / (2 * cpp))) ** 2


def V_int(r, theta, t):
    """
    Laser–atom interaction potential in the length gauge.

    :param r: Radial coordinate in atomic units (float).
    :param theta: Polar angle in radians (float).
    :param t: Time in atomic units (float).
    :return: Interaction potential at (r, θ, t) (float).

    Reference
    ----------
    - Section 2.3 -- *Time Evolution of atomic wavefunction interacting with external strong-field LASER
    """
    return -E_field(t) * r * np.cos(theta)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    Absorbing function                    |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Absorber_func(r):
    """
    Radial absorber function.

    :param r: Radial coordinate (a.u.).
    :return: Absorber value at r (float).

    Reference
    ----------
    - Section 2.3.9 -- *The Absorber mask function and Absorbing layer*
    """
    if 0 < r <= r0: return 1
    elif r0 < r < r_max:
        return np.cos(np.pi * (r - r0) / (2 * (r_max - r0))) ** 0.25



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      Single-active-electron (SAE) model potentials       |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def potential_V_SAE_M1(r, atom='Ne'):
    """
    Single-active-electron (SAE) model potential (Model-1).

    :param r: Radial coordinate (a.u.).
    :param atom: Atomic species label (default: 'Ne').
    :return: SAE model potential V(r) (a.u.).
    :raises ValueError: If atom is not found in the parameter table.

    Reference
    ----------
    - Section 3.1.1 -- *Atomic model potential: SAE-M1*
    """
    params = atomic_params_SAE_M1.get(atom)
    if params is None: raise ValueError(f"Atom or ion '{atom}' not found in table.")
    Zc, a1, a2, a3, a4, a5, a6 = params["Zc"], params["a1"], params["a2"], params["a3"], params["a4"], params["a5"], params["a6"]
    return -(Zc + a1*np.exp(-a2*r) + a3*r*np.exp(-a4*r) + a5*np.exp(-a6*r)) / r

def potential_V_SAE_M2(r, atom='Ne'):
    """
    Single-active-electron (SAE) model potential (Model-2).

    :param r: Radial coordinate (a.u.).
    :param atom: Atomic species label (default: 'Ne').
    :return: SAE model potential V(r) (a.u.).
    :raises ValueError: If atom is not found in the parameter table.

    Reference
    ----------
    - Section 3.1.2 -- *Atomic model potential: SAE-M2*
    """
    params = atomic_params_SAE_M2.get(atom)
    if params is None: raise ValueError(f"Atom or ion '{atom}' not found in table.")

    C0, Zc, c = params['C0'], params['Zc'], params['c']
    a1, a2, a3 = params['a1'], params['a2'], params['a3']
    b1, b2, b3 = params['b1'], params['b2'], params['b3']

    V_long = -C0 / r
    V_short = -Zc * np.exp(-c * r) / r
    V_shell = (a1 * np.exp(-b1 * r)) + (a2 * np.exp(-b2 * r)) + (a3 * np.exp(-b3 * r))
    return V_long + V_short - V_shell       # overall -ve sign in V_shell.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Helpful functions that simplify life           |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def V_eff(l, x):
    """
    Effective potential for a particle in a central field.

    :param l: Angular momentum quantum number (int).
    :param x: Radial coordinate (float).
    :return: Effective potential value at x (float).
    """
    return -1 / x + l*(l+1) / (2*x**2)


def P_N(x):
    """
    Legendre polynomial of degree N at x.

    :param x: Point(s) of evaluation (float or array-like).
    :param N: Polynomial degree (int).
    :return: Value(s) of P_N(x) (float or numpy.ndarray).
    """
    return legendre(N)(x)


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


def generate_states(l):
    """
    Generate a list of electronic states for a given orbital angular momentum quantum number.

    The states are labeled by the principal quantum number (starting from :math:`n = \ell + 1`)
    followed by the spectroscopic orbital letter (s, p, d, f, g, ...).

    Example:

    >>> generate_states(1)
    ['2p', '3p', '4p', ..., '199p']

    :param l: Orbital angular momentum quantum number :math:`\ell` (int).
              Supported up to :math:`\ell = 10` (m orbital).
    :return: List of state labels as strings (list of str).
    """
    orbital_types = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h',
                     6: 'i', 7: 'j', 8: 'k', 9: 'l', 10: 'm'}
    return [str(i) + orbital_types.get(l, '') for i in range(l + 1, 200)]


def state_name(n, l):
    """
    Converts quantum numbers n and l to a string representation of the atomic state.

    :param n: Principal quantum number (n ≥ 1).
    :param l: Orbital angular momentum quantum number (l ≥ 0).
    :return: Atomic state string (e.g., '1s', '2p', '3d').
    :raises ValueError: If l is outside the allowed range (0–6).

    Example:

    >>> state_name(1, 0)
    '1s'
    """
    orbital_letters = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i'}
    if l in orbital_letters:
        return f"{n}{orbital_letters[l]}"
    else:
        raise ValueError(f"Invalid orbital angular momentum quantum number l={l}. Allowed values are 0 to 6.")




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 Printing info                  |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ != '__main__':
    print('~~~~~~~~~~~~~: Grid info :~~~~~~~~~~~~~')
    print('mapping param (L_map)       :', L_map)
    print('mapping param (r_max)       :', r_max)
    print('radial colloc points (N-1)  :', N-1)
    print('angular colloc points (L+1) :', L+1, '\n')

    # print('~~~~~~~~~~~~~~: Spectra :~~~~~~~~~~~~~~')  # TODO: move to evolution.py
    # print(f'Ip (a.u)              : {Ip:.3f}')
    # print(f'Up (a.u)              : {Up_au:.5f}')
    # print(f'N_cutoff              : {N_cut:.3f}')
    # print(f'Keldysh parameter (γ) : {Keldysh(Ip, Up_au):.3f}\n')

    print('~~~~~~~~~: Atom & Laser info :~~~~~~~~~')
    print(f'atom system   : {evolving_atom}')
    print(f'initial state : {state_name(n, l)}')
    print(f'I0 (W/cm2)    : {I0:.2e}')
    print('I0 (a.u)      :', I0 / Int_0)
    print('E0 (a.u)      :', E0_au)
    print('λ  (nm)       :', lambda_nm)
    print('λ  (a.u)      :', lambda_nm * 10 ** -9 / a0)
    print('w0 (a.u)      :', w0)
    print('T (a.u)       :', T)
    print('T (f.s)       :', T * T0 * 10**15)
    print('tf (a.u)      :', tf)
    print('tf (f.s)      :', tf * T0 * 10**15)
    print('dt (a.u)      :', dt)
    print('dt (atto)     :', dt * T0 * 10**18)
    print('nopt          :', len(t), '\n')

    print('~~~~~~~~~~~: S-matrix info :~~~~~~~~~~~')
    print('l_max         :', l_max)
    print('k_max         :', k_max)
    print('dt            :', dt, '\n')