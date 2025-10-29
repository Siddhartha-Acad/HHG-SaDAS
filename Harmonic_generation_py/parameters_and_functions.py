"""
File: parameters_and_functions.py
Project: HHG-SaDAS

Code Description:
    - Contains all parameters that define the entire system and numerical requirements.
    - All other simulation codes fetch parameters and functions from this script.

*** [NAMING CONVENTION FOLLOWED IN THE ENTIRE PACKAGE] ***

Quantum Number Convention:
--------------------------
- n here is the radial index for a given l.
- For a given orbital angular momentum quantum number l ≥ 0,
  the **actual principal quantum number** is: n_effective = n + l
- This mapping corresponds to standard spectroscopic notation.

Example mapping of quantum numbers (n, l, m) to standard orbitals:
------------------------------------------
  l   |   n   |   m   |   Orbital / State
  ----------------------------------------
  0   |   1   |   0   |    1s
  1   |   1   |   0   |    2p_{z}
  1   |   2   |   0   |    3p_{z}
  1   |   1   |   1   |    2p_{x}
  1   |   2   |   1   |    3p_{x}
  3   |   1   |   0   |    4f_{z^3}
------------------------------------------

Key Ideas:
-----------
1. n is the **radial index**, starting from 1 for each l.
2. l determines the **orbital type**: 0→s, 1→p, 2→d, 3→f, etc.
3. m determines the **magnetic sublevel**, i.e., orbital orientation.
4. The **actual principal quantum number** = n + l.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

Notes:
------
- This file is part of the HHG-SaDAS package, developed during the MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
"""

import warnings
from pathlib import Path
from scipy.special import legendre
from scipy.special import factorial, lpmv
from Atomic_units import Int_0, omega_au, T0
from Harmonic_generation_py.conf_model_bank import *

this_dir = Path(__file__).resolve().parent     # Relative path system
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
print('*** RuntimeWarning     : Blocked from parameters_and_functions.py ***')
print('*** DeprecationWarning : Blocked from parameters_and_functions.py ***\n')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Atom, SAE and Confinement            |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n = 1; l = 0; m = 0        # defines initial state. [NOTE]: PLEASE ADAPT TO THE NAMING CONVENTION MENTIONED IN TEH DOCSTRING.
evolving_atom = 'H'        # Atoms are listed down in 'SAE dataset' section.
SAE_model = 'SAE-M1'       # Single active electron model; option: SAE_model = 'SAE-M1' or 'SAE-M2'. [NOTE]: For 'Xe' always use 'SAE-M1'

confined = False                    # whether the atom is confined or not?
confinement_model = 'P-Gau'         # which type of confinement potential? Options: 'ASW', 'GASW', 'Lor', 'SSW', 'Gau', 'P-Gau'
save_Egvals_with_Smatrix = True     # Eigenvalues for all l=(m, l_max+m) will be saved in 'GPSM_states_S-matrix/GPSM_states_and_Smatrix_data/
                                    # setting: save_Egvals_with_Smatrix=True will save eigenvalues that are used to make an energy level diagram
                                    # in: HHG-SaDAS/Harmonic_generation_py/GPSM_states_S-matrix/Energy_level_diagram.py


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                GPSM Parameters                 |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = 200                     # P'_N(xj) = 0 ; radial grid size: len(colloc_pt) = N-1
L = 20                      # must be >= l ; angular grid size: len(theta_k) = L+1
l_max = 20                  # Number of partial waves = number of S-matrices = l_max+1
k_max = 50                  # number of GPSM states (maximum k index) in S matrix
L_map = 80; r_max = 200     # radial mapping parameters
r0 = 150                    # absorber layer thickness = (r_max - r0) a.u.


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
#   Time evolution controls: time_evolution.py   |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eta_t = 0.027               # Execution time for a single time-step (dt) evolution. eta_t = 0.03 is the execution speed achieved on my system.
time_step = 1000            # number of time steps desired for evolution. Maximum possible steps = len(t)-1. {-1 because time_step used as index}
show_E_field = False        # Whether to display the laser electric field before the evolution starts. (plot will remain open until you kill it).
print_serial_prog = True    # when True, running time_evolution.py will print progress. Example:  {Evolution step 49    : 50.00%}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      collocation grid (radial & angular)       |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colloc_file = f'Algo-3_N={N}_Gauss_Lobatto_collocation_points.txt'                      # File that holds the collocation points.
colloc_file = this_dir / 'collocation_points' / 'generator' / colloc_file               # fetching collocation data from relative path.
colloc_pt = np.loadtxt(colloc_file, skiprows=1, usecols=0)                              # Gauss-Lobatto collocation points.

int_w = 2 / (N * (N + 1) * (legendre(N)(colloc_pt))**2)           # Gauss-Lobatto  Quadrature weights: w_j
roots, weights = np.polynomial.legendre.leggauss(L+1)             # Gauss-Legendre Quadrature weights and collocation points (or, nodes): x_k
theta_k = np.arccos(roots)                                        # Angular collocation points: cos(theta_k) = x_k


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
    alpha_map = 2 * Lmap / r_max
    return Lmap * (1 + x) / (1 - x + alpha_map)


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
    alpha_map = 2 * Lmap / r_max
    return Lmap * (alpha_map + 2) / (1 - x + alpha_map)**2



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
def H(l_val, i, j, model=SAE_model):
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

    :param l_val: Angular momentum quantum number :math:`\ell` (int).
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

    conf_term = conf_selector(confinement_model, f(colloc_pt[i]))[0] if confined else 0.0

    if model == 'SAE-M1':
        term2 = (l_val * (l_val + 1) / (2 * f(colloc_pt[i]) ** 2) +
                 potential_V_SAE_M1(f(colloc_pt[i]), atom=evolving_atom) +
                 conf_term)
        return term1 + term2
    elif model == 'SAE-M2':
        term2 = (l_val * (l_val + 1) / (2 * f(colloc_pt[i]) ** 2) +
                 potential_V_SAE_M2(f(colloc_pt[i]), atom=evolving_atom) +
                 conf_term)
        return term1 + term2



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Constant factors -- associated leg. poly -- Spherical Harmonics        |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def N_fact(l_val, m_val):
    """
    Normalization constant of Y_lm
    """
    return (-1)**m_val * np.sqrt((2 * l_val + 1) * factorial(l_val - m_val) / (4 * np.pi * factorial(l_val + m_val)))


def C_fact(l_val, m_val):
    """
    Orthogonality constant factor of P_lm
    """
    return 2 * factorial(l_val + m_val) / ((2 * l_val + 1) * factorial(l_val - m_val))


def a_legendre(l_val, m_val, x):
    """
    Associated Legendre polynomial without Condon-Shortley phase
    """
    return lpmv(m_val, l_val, x) * (-1)**m_val


def Y_lm(l_val, m_val, x):
    """
    Y_lm(x) = N_lm * P_lm(x)
    """
    return N_fact(l_val, m_val) * a_legendre(l_val, m_val, x)


def Y_lm_array(l_max, m, roots):
    """
    Fully vectorized computation of Y_lm(l, m, x) for all l in [0, l_max] and x in roots.
    """
    l_vals = np.arange(0, l_max + 1)[:, None]    # shape (l_max+1, 1)
    x = np.asarray(roots)[None, :]               # shape (1, L+1)

    # Broadcast l_vals and x to same shape
    # lpmv can handle broadcasting since m is scalar
    P_lm = lpmv(m, l_vals, x)                    # shape (l_max+1, L+1)

    # Normalization factors
    N_lm = (-1)**m * np.sqrt(
        (2 * l_vals + 1) *
        factorial(l_vals - m, exact=False) /
        (4 * np.pi * factorial(l_vals + m, exact=False))
    )

    return N_lm * P_lm



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  Partial-wave g_lm(r) calculating function                   |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def g_lm(Psi_t, glm_arr):
    r"""
    Compute radial partial-wave projections :math:`g_{\ell m}(r)` from an angular–radial
    wavefunction using Gauss–Legendre quadrature, with optional memory reuse.

    .. math::
        g_{\ell m}(r) = \frac{1}{N_{\ell m} \cdot C_{\ell m}} \sum_k w_k \, P_\ell^{(m)}(\cos \theta_k) \, \psi(\theta_k, r)

    where :math:`w_k` are Gauss–Legendre quadrature weights,
    :math:`P_\ell^{(m)}(\cos \theta_k)` are precomputed associated Legendre polynomials,
    :math:`N_{\ell m}` and :math:`C_{\ell m}` are normalization factors from ``N_fact`` and ``C_fact``.

    :param Psi_t: Angular–radial wavefunction :math:`\psi(\theta, r)` at a fixed time,
                  evaluated at Gauss–Legendre quadrature points.
                  Shape = ``(n_theta, N-1)``, where
                  ``Psi_t[k, i] = ψ(θ_k, r_i)``.
    :type Psi_t: ndarray

    :param glm_arr: Optional output array for storing results in-place to avoid
                   repeated allocations. If None, a new array is created.
                   Shape = ``(l_max+1, N-1)``.
    :type glm_arr: ndarray or None

    :return: Radial partial-wave projections :math:`g_{\ell m}(r)` for
             :math:`\ell = 0, \ldots, \ell_\text{max}` at fixed magnetic quantum number :math:`m`.
             Each row corresponds to one angular momentum component :math:`\ell`.
    :rtype: ndarray, shape ``(l_max+1, N-1)``

    Notes
    -----
    This function computes the partial-wave expansion coefficients for a specific
    magnetic quantum number :math:`m`. The associated Legendre polynomials
    :math:`P_\ell^{(m)}(\cos \theta)` are precomputed in ``a_legendre_vals`` using
    the ``a_legendre`` function.

    This function is optimized for performance-critical use cases
    (e.g. propagation loops), where avoiding memory reallocation
    at each step is important.

    Example
    -------
    >> glm_empty = np.empty((l_max + 1, N-1), dtype=np.complex128)
    >> for ti in range(len(t)-1):
    ...     # Psi_t = updated wavefunction
    ...     gl_vals = g_lm(Psi_t, glm_empty)  # memory reused each iteration

    References
    ----------
    - Appendix D -- *Derivation of the general partial wave formula*
    - Section 2.3.4 -- *The Partial wave expansion*
    - Section 3.3.2 -- *The partial-wave expansion problem*
    - Section 3.3.3 -- *Generalization over magnetic quantum number m*
    """
    for l_ind in range(l_max+1):
        glm_arr[l_ind] = np.tensordot(weights * a_legendre_vals[l_ind], Psi_t, axes=([0], [0])) / (N_fact(l_ind+m, m) * C_fact(l_ind+m, m))
    return glm_arr

a_legendre_vals = np.array([[a_legendre(l_index+m, m, root) for root in roots] for l_index in range(l_max+1)])  # shape = (l_max+1, L+1)
vect_norm_fact = np.array([1.0 / (N_fact(l_index+m, m) * C_fact(l_index+m, m)) for l_index in range(l_max+1)])  # shape = (l_max+1, )
weighted_legendre_vals = weights[None, :] * a_legendre_vals                                                     # shape = (l_max+1, L+1)

def g_lm_vect(Psi_t, glm_arr):
    r"""
    Compute the radial partial-wave projections :math:`g_{\ell m}(r)` in a fully
    vectorized manner using precomputed Legendre polynomials and normalization
    constants.

    This function is mathematically equivalent to ``g_lm()``, but it eliminates
    all Python-level loops and redundant intermediate allocations for maximum
    performance. It performs the projection over all angular momentum quantum
    numbers :math:`\ell` in a single matrix multiplication using BLAS-optimized
    NumPy routines.

    The underlying formula is identical to the standard partial-wave expansion:

    .. math::
        g_{\ell m}(r) =
            \frac{1}{N_{\ell m} C_{\ell m}}
            \sum_k w_k P_\ell^{(m)}(\cos \theta_k) \, \psi(\theta_k, r)

    where:
        - :math:`w_k` are Gauss–Legendre quadrature weights,
        - :math:`P_\ell^{(m)}(\cos \theta_k)` are precomputed associated Legendre polynomials,
        - :math:`N_{\ell m}` and :math:`C_{\ell m}` are normalization factors.

    Parameters
    ----------
    Psi_t : ndarray, shape (n_theta, N-1)
        Angular–radial wavefunction ψ(θ, r) at a fixed time step,
        evaluated at Gauss–Legendre quadrature points.

    glm_arr : ndarray, shape (l_max+1, N-1)
        Preallocated output array for in-place storage of
        partial-wave projections g_{ℓm}(r). This avoids repeated
        memory allocations in time propagation loops.

    Returns
    -------
    glm_arr : ndarray, shape (l_max+1, N-1)
        Radial partial-wave projections for ℓ = 0, …, ℓ_max
        corresponding to the given magnetic quantum number m.

    Notes
    -----
    - This version is **fully vectorized** and performs all ℓ projections
      simultaneously via a single matrix multiplication.
    - Functionally identical to :func:`g_lm`, differing only in implementation.
    - Significantly faster and more memory-efficient for long time evolutions.

    Example
    -------
    >>> glm_arr = np.empty((l_max + 1, N - 1), dtype=np.complex128)
    >>> g_lm_vect(Psi_t, glm_arr)  # In-place computation
    """
    np.dot(weighted_legendre_vals, Psi_t, out=glm_arr)
    glm_arr *= vect_norm_fact[:, None]
    return glm_arr



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  Partial-wave evolution                  |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def G(Sl_matrix, gl_arr, l_ind):
    r"""
    Compute the third bracket of Eq.(2.86),
    which is a matrix multiplication of the S-matrix block with the partial-wave projection.

    .. math::
        G_\ell(r_i, t) = \sum_{j=1}^{N-1} S_{ij}(\ell) \, g_\ell(r_j, t)

    Parameters
    ----------
    Sl_matrix : ndarray, shape (l_max+1, N-1, N-1)
        S-matrix blocks :math:`S_\ell` for each angular momentum :math:`\ell`.
    gl_arr : ndarray, shape (l_max+1, N-1)
        Radial partial-wave projections :math:`g_\ell(r, t)`.
    l_ind : int
        Angular momentum index :math:`\ell`.

    Returns
    -------
    ndarray, shape (N-1,)
        Transformed radial component :math:`G_\ell(r, t)`.
    """
    return np.dot(Sl_matrix[l_ind], gl_arr[l_ind])



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      dipole moment                       |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def alpha_lm(l_val, m_val):
    r"""
    This factor appears in the expression of general dipole moment formula (see Eq.~E.18)

    .. math::
        \alpha_{l,m} = \sqrt{\frac{(l + 1)^2 - m^2}{(2l + 1)(2l + 3)}}

    Parameters
    ----------
    l_val : int or ndarray
        Angular momentum quantum number(s) :math:`\ell`. Can be a scalar or a NumPy array.
    m_val : int or ndarray
        Magnetic quantum number(s) :math:`m` (:math:`|m| \leq \ell`). Can be a scalar or a NumPy array.
        If one is an array and the other is a scalar, broadcasting is applied.

    Returns
    -------
    float
        The computed dipole factor :math:`\alpha_{l,m}`.
    """
    return np.sqrt(((l_val + 1) ** 2 - m_val ** 2) / ((2 * l_val + 1) * (2 * l_val + 3)))


def dipole_moment(r, glm_arr):
    r"""
    Computes the dipole moment for a given radial grid and a set of partial-waves at an instant of time.

    .. math::
        d_{i\sigma}(t) = 2 \sum_\ell \alpha_{\ell,m} \int r \, \mathrm{Re}\left[g_\ell^*(r, t) g_{\ell+1}(r, t)\right] \, dr

    Parameters
    ----------
    r : ndarray, shape (N-1,)
        Nonlinear Radial collocation grid points: r(x)
    glm_arr : ndarray, shape (l_max+1, N-1)
        Partial-wave amplitudes :math:`g_\ell(r)` of the wavefunction.

        - `l_max+1` : number of angular momentum channels.
        - `N-1`     : number of radial grid points.

    Returns
    -------
    float
        Total dipole moment computed from the partial-wave amplitudes.

    Notes
    -----
    - `l_ind_arr = np.arange(l_max)` is used to iterate over angular momentum indices : precomputed
    - `alpha_factor = alpha(l_ind_arr + m, m)` : precomputed
    """
    integrals = np.array([np.sum(r * np.real(np.conj(glm_arr[l_ind]) * glm_arr[l_ind + 1])) for l_ind in l_ind_arr])
    return 2 * np.sum(alpha_factor * integrals)

l_ind_arr = np.arange(l_max)            # l_max because in Eq.~E.21, the `l' index goes from (m) to (l_max+m-1). So, in total (l_max-1).
alpha_factor = alpha_lm(l_ind_arr+m, m)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   survival probability                   |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Ps(glm_arr):
    r"""
    Calculates the total survival probability at a given instant.

    Parameters
    ----------
    glm_arr : ndarray, shape (l_max+1, N-1)
        Array of radial partial-wave amplitudes :math:`g_\ell(r)`.

        - `l_max+1` : number of angular momentum channels.
        - `N-1`     : number of radial grid points.

    Returns
    -------
    float
        The survival probability :math:`P_s = \sum_\ell \sum_r |g_\ell(r)|^2`.

    Reference
    ----------
    - Section 2.4.3 -- *Correlation function and Survival probability*
    """
    return np.sum(np.abs(glm_arr)**2)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ponderomotive force;  Keldysh parameter; Harmonic cut-off |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Up(field_amp_au, omega0_au):
    r"""
    Ponderomotive energy :math:`U_p`.

    .. math::
        U_p = \frac{E_0^2}{4 \omega_0^2}

    :param field_amp_au: Electric field amplitude (a.u.).
    :param omega0_au: Angular frequency of the laser field (a.u.).
    :return: Ponderomotive energy (a.u.).
    """
    return field_amp_au**2 / (4 * omega0_au**2)


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


def N_cutoff(Ip_au, Up_au, omega0_au):
    r"""
    Harmonic cut-off order :math:`N_c`.

    .. math::
        N_c = \frac{I_p + 3.17 U_p}{\omega_0}

    :param Ip_au: Ionization potential (a.u.).
    :param Up_au: Ponderomotive energy (a.u.).
    :param omega0_au: Angular frequency of the laser field (a.u.).
    :return: Cut-off harmonic order (dimensionless).
    """
    return (Ip_au + 3.17 * Up_au) / omega0_au



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      LASER electric field and interaction-potential      |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def E_field(t_val):
    r"""
    Laser electric field with a sine-squared envelope.

    .. math::
        E(t) = E_0 \, \sin(\omega_0 t) \,
        \left[ \sin\!\left(\frac{\omega_0 t}{2 N_c}\right) \right]^2

    :param t_val: Time (a.u.).
    :return: Electric field amplitude at time t (a.u.).

    References
    ----------
    - Section 2.3 — *Time evolution of the atomic wavefunction interacting with an external strong-field laser*
    """
    return E0_au * np.sin(w0 * t_val) * (np.sin(w0 * t_val / (2 * cpp))) ** 2


def V_int(r, theta, t_val):
    """
    Laser–atom interaction potential in the length gauge.

    :param r: Radial coordinate in atomic units (float).
    :param theta: Polar angle in radians (float).
    :param t_val: Time in atomic units (float).
    :return: Interaction potential at (r, θ, t) (float).

    Reference
    ----------
    - Section 2.3 -- *Time Evolution of atomic wavefunction interacting with external strong-field LASER
    """
    return -E_field(t_val) * r * np.cos(theta)



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
#              Confinement potential selector              |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def conf_selector(input_model, r):
    """
    :param input_model: String specifying the confinement model.
                        Options: 'ASW', 'GASW', 'Lor', 'SSW', 'Gau', 'P-Gau'
    :param r: Radial coordinate (float or array-like).
    :return: Tuple (v_pot, info_string)
             - v_pot: Potential values as a NumPy array.
             - info_string: String with parameter values.
    """
    if input_model == 'ASW':
        return V_ASW(r)
    elif input_model == "GASW":
        return V_GASW(r)
    elif input_model == 'Lor':
        return V_Lorentz(r)
    elif input_model == 'SSW':
        return V_SSW(r)
    elif input_model == 'Gau':
        return V_Gaussian(r)
    elif input_model == 'P-Gau':
        return V_PowerExpo(r)
    else:
        raise ValueError(f"Unknown confinement model: {input_model}")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Helpful functions that simplify life           |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def V_eff(l_val, x):
    """
    Effective potential for a particle in a central field.

    :param l_val: Angular momentum quantum number (int).
    :param x: Radial coordinate (float).
    :return: Effective potential value at x (float).
    """
    return -1 / x + l_val*(l_val + 1) / (2 * x ** 2)


def P_N(x):
    """
    Legendre polynomial of degree N at x.

    :param x: Point(s) of evaluation (float or array-like).
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


def generate_states(l_val):
    r"""
    Generate a list of electronic states for a given orbital angular momentum quantum number.

    The states are labeled by the principal quantum number (starting from :math:`n = \ell + 1`)
    followed by the spectroscopic orbital letter (s, p, d, f, g, ...).

    Example:

    >> generate_states(1)
    ['2p', '3p', '4p', ..., '199p']

    :param l_val: Orbital angular momentum quantum number :math:`\ell` (int).
              Supported up to :math:`\ell = 10` (m orbital).
    :return: List of state labels as strings (list of str).
    """
    orbital_types = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h',
                     6: 'i', 7: 'j', 8: 'k', 9: 'l', 10: 'm'}
    return [str(i) + orbital_types.get(l_val, '') for i in range(l_val + 1, 200)]


def state_name(n_val, l_val):
    r"""
    Converts quantum numbers n and l to a string representation of the atomic state.

    :param n_val: Principal quantum number (n ≥ 1).
    :param l_val: Orbital angular momentum quantum number (l ≥ 0).
    :return: Atomic state string (e.g., '1s', '2p', '3d').
    :raises ValueError: If l is outside the allowed range (0–6).

    Example:

    >>> state_name(1, 0)
    '1s'
    """
    orbital_letters = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i'}
    if l_val in orbital_letters:
        return f"{n_val}{orbital_letters[l_val]}"
    else:
        raise ValueError(f"Invalid orbital angular momentum quantum number l={l_val}. Allowed values are 0 to 6.")




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 Printing info                  |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ != '__main__':
    print('~~~~~~~~~~~~~: Grid info :~~~~~~~~~~~~~')
    print('mapping param (L_map)       :', L_map)
    print('mapping param (r_max)       :', r_max)
    print('radial colloc points (N-1)  :', N-1)
    print('angular colloc points (L+1) :', L+1, '\n')

    # print('~~~~~~~~~~~~~~: Spectra :~~~~~~~~~~~~~~')
    # print(f'Ip (a.u)              : {Ip:.3f}')
    # print(f'Up (a.u)              : {Up_au:.5f}')
    # print(f'N_cutoff              : {N_cut:.3f}')
    # print(f'Keldysh parameter (γ) : {Keldysh(Ip, Up_au):.3f}\n')

    print('~~~~~~~~~: Atom & Laser info :~~~~~~~~~')
    if not confined:
        print(f'atom system   : {evolving_atom}')
    else:
        print(f'atom system   : {evolving_atom}@C60')
        print(f'conf. model   : {confinement_model}')

    print(f'initial state : ({n=}, {l=}, {m=}) ~ {state_name(n+l, l)} --> time_evolution.py')      # PRINCIPLE QUANTUM NUMBER = n+l
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
