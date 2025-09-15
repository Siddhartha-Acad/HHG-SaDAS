"""
File: Simplified_alpha_lm
Project: HHG-SaDAS
Code Description:
    | [just for demo]
    | Computes and verifies the prefactor α_lm, which appears in the generalized
    | dipole moment formula (see Eq. E.18 in my thesis).
    | The original expression for alpha_lm was quite complicated; this code tests
    | the correctness of the simplified version of alpha_lm. The supporting calculations
    | are also given in this directory.

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

import math

def factorial(n):
    return math.factorial(n)

def N_lm(l, m):
    num = (2*l + 1) * factorial(l - m)
    den = 4 * math.pi * factorial(l + m)
    return (-1)**m * math.sqrt(num / den)

def C_lm(l, m):
    num = 2 * factorial(l + m)
    den = (2*l + 1) * factorial(l - m)
    return num / den

def alpha_lm_original(l, m):
    """From original long expression"""
    term1 = 2 * math.pi * (l + m + 1) / (2*l + 3)
    Clm = C_lm(l, m)
    Nlm = N_lm(l, m)
    Nlm1 = N_lm(l + 1, m)
    return term1 * Clm * Nlm * Nlm1

def alpha_lm_simplified(l, m):
    """From final simplified expression"""
    num = (l + 1)**2 - m**2
    den = (2*l + 1) * (2*l + 3)
    return math.sqrt(num / den)

# Test for some values of l and m
print(f"{'l':<3} {'m':<3} {'Original alpha_lm':<20} {'Simplified alpha_lm':<20} {'Difference'}")
for l in range(1, 6):
    for m in range(0, l + 1):  # valid m: 0 ≤ m ≤ l
        orig = alpha_lm_original(l, m)
        simp = alpha_lm_simplified(l, m)
        diff = abs(orig - simp)
        print(f"{l:<3} {m:<3} {orig:<20.12f} {simp:<20.12f} {diff:.2e}")
