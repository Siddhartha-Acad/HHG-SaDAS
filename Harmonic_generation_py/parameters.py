"""
File: parameters.py
Project: HHG-SaDAS

Code Description:
    - Contains all parameters that define the entire system and numerical requirements.
    - All other simulation codes fetch parameters from this script.

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
import numpy as np
from pathlib import Path
from scipy.special import legendre
from scipy.special import factorial, lpmv
from Atomic_units import Int_0, omega_au, T0

this_dir = Path(__file__).resolve().parent     # Relative path system
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


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
total_states = 5            # how many states you want to keep in the GPSM_state file (.dat)

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
eta_t = 0.00154             # Execution time for a single time-step (dt) evolution. eta_t = 0.00154 is the execution speed achieved on my system.
time_step = len(t) - 1      # number of time steps desired for evolution. Maximum possible steps = len(t)-1. {-1 because time_step used as index}
show_E_field = False        # Whether to display the laser electric field before the evolution starts. (plot will remain open until you kill it).
print_serial_prog = True    # when True, running time_evolution.py will print progress. Example:  {Evolution step 49    : 50.00%}


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

