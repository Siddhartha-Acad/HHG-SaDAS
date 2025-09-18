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
