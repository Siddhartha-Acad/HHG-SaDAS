"""
File: parameters.py
Project: HHG-SaDAS
Code Description:
    This contains all parameters and functions that defines the entire system and numerical requirements.
    All other simulating codes, fetch parameters from this script.


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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                 GPSM Parameters                |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N = 200                     # P'_N(xj) = 0 ; radial grid size: len(colloc_pt) = N-1
L = 20                      # must be >= l ; number of S matrix gl(r) in partial wave expansion.
k_max = 50                  # number of GPSM states (maximum k index) in S matrix
L_map = 20 ; r_max = 200    # radial mapping parameters
r0 = 150                    # absorber layer thickness: (r_max - r0) a.u.