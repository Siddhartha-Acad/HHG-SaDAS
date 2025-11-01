"""
File: conf_model_bank.py
Project: HHG-SaDAS
Code Description:
    | Library of Confinement Potential Models used for simulating
    | harmonic generation in confined noble-gas atoms.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Contains implementations of confinement potentials:
  * Annular Square Well (ASW) potential
  * Gaussian (Gau) potential
  * Gaussian Annular Square Well (GASW) potential
  * Lorentzian (Lor) potential
  * Smooth Square Well (SSW) potential
  * Power-Gaussian (P-Gau) potential

- Each function returns the potential values and a parameter info string.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import numpy as np
from Atomic_units import a0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Annular Square Well (ASW) potential            |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def V_ASW(r):
    """
    :param r: Radial coordinate (float or array-like).
    :return: Tuple (v_asw, info_string)
             - v_asw: Potential values as a NumPy array.
             - info_string: String with parameter values.
    """
    ASW_U_au = 0.56
    ASW_Delta_au = 2.8
    ASW_rc_au = 6.7

    r = np.atleast_1d(r)
    mask = (ASW_rc_au - ASW_Delta_au / 2 <= r) & (r <= ASW_rc_au + ASW_Delta_au / 2)
    v_asw = np.zeros_like(r, dtype=float)
    v_asw[mask] = -ASW_U_au

    ASW_info_string = f'ASW_U={ASW_U_au}_Delta={ASW_Delta_au}_rc={ASW_rc_au}'
    return v_asw, ASW_info_string


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 Gaussian (Gau) potential                 |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def V_Gaussian(r):
    """
    :param r: Radial coordinate (float or array-like).
    :return: Tuple (v_gau, info_string)
             - v_gau: Gaussian potential values as a NumPy array.
             - info_string: String with parameter values.
    """
    w0_Ryd = 0.647
    sigma_A = 0.57
    rc_A = 3.54

    w0_au = 0.5 * w0_Ryd; sigma_au = sigma_A * 10 ** -10 / a0; rc_au = rc_A * 10 ** -10 / a0
    Gau_info_string = f'Gau_w0={w0_au:.3f}_sigma={sigma_au:.3f}_rc={rc_au:.3f}'
    return -w0_au * np.exp(-(r - rc_au)**2 / sigma_au**2), Gau_info_string


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      Gaussian Annular Square Well (GASW) potential       |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def V_GASW(r):
    """
    :param r: Radial coordinate (float or array-like).
    :return: Tuple (v_gasw, info_string)
             - v_gasw: Combined Gaussian + ASW potential values as a NumPy array.
             - info_string: String with parameter values.
    """
    p_rat = 1.2        # p in thesis. using p_rat (p ratio) to avoid confusion with P-Gau
    depth_w0 = 0.324
    delta_au = 2.8
    sigma_au = 1.70
    rc_au = 6.69

    amp_Gau = depth_w0 / (1 + 1 / p_rat)
    amp_ASW = amp_Gau / p_rat

    r = np.atleast_1d(r)
    V_Gau = -amp_Gau * np.exp(-((r - rc_au) / (np.sqrt(2) * sigma_au))**2)
    mask = (rc_au - delta_au / 2 <= r) & (r <= rc_au + delta_au / 2)
    V_ASW2 = np.zeros_like(r, dtype=float)
    V_ASW2[mask] = -amp_ASW
    v_gasw = V_Gau + V_ASW2

    info_string = f'GASW_prat={p_rat}_w0={depth_w0:.3f}_sigma={sigma_au:.2f}_delta={delta_au}_rc={rc_au:.2f}'
    return v_gasw, info_string


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                Lorentzian (Lor) potential                |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def V_Lorentz(r):
    """
    :param r: Radial coordinate (float or array-like).
    :return: Tuple (v_lorentz, info_string)
             - v_lorentz: Lorentzian potential values as a NumPy array.
             - info_string: String with parameter values.
    """
    U_au = 0.04663
    A_au = -1.11903
    W_au = 1.86046
    rc_au = 6.68963

    info_string = f'Lor_U={U_au:.3f}_A={A_au:.3f}_W={W_au:.3f}_rc={rc_au:.3f}'
    return U_au + 2 * A_au * W_au / (np.pi * (4 * (r - rc_au)**2 + W_au**2)), info_string


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            Smooth Square Well (SSW) potential            |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def V_SSW(r):
    """
    :param r: Radial coordinate (float or array-like).
    :return: Tuple (v_ssw, info_string)
             - v_ssw: Smooth square-well potential values as a NumPy array.
             - info_string: String with parameter values.
    """
    w0_au = 0.322362
    delta_au = 2.065270
    rc_au = 6.689630
    k_smooth = 14.288530   # Sharpness of transition

    def sigmoid(x, k):
        return 1 / (1 + np.exp(-k * x))

    # Smooth square-well potential
    S1 = sigmoid(r - (rc_au - delta_au / 2), k_smooth)
    S2 = sigmoid(r - (rc_au + delta_au / 2), k_smooth)

    SSW_info_string = f'SSW_w0={w0_au:.3f}_delta={delta_au:.3f}_k={k_smooth:.3f}_rc={rc_au:.3f}'
    return -w0_au * (S1 - S2), SSW_info_string


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#             Power-Gaussian (P-Gau) potential             |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def V_PowerExpo(r):
    """
    :param r: Radial coordinate (float or array-like).
    :return: Tuple (v_pexp, info_string)
             - v_pexp: Power-exponential potential values as a NumPy array.
             - info_string: String with parameter values.
    """
    p = 10
    w0_Ryd = 0.647
    sigma_A = 0.57
    rc_A = 3.54

    w0_au = 0.5 * w0_Ryd; sigma_au = sigma_A * 10 ** -10 / a0; rc_au = rc_A * 10 ** -10 / a0
    PGau_info_string = f'Pexp_p={p}_w0={w0_au:.3f}_sigma={sigma_au:.3f}_rc={rc_au:.3f}'
    return -w0_au * np.exp(-(r - rc_au)**p / sigma_au**p), PGau_info_string

