"""
File: keldysh_param_SI_directly.py
Project: HHG-SaDAS
Code Description:
    | Computes the Keldysh parameter for an atom under a laser field.
    | The calculation uses fundamental constants, laser parameters, and the
    | ionization potential to evaluate whether the interaction is in the
    | tunneling or multiphoton regime.

Author: Siddhartha Mithiya
Affiliation: IIT Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Input laser wavelength: 1064 nm.
- Input intensity: 5e13 W/cm^2.
- Ionization potential used: Hydrogen 1s, 13.6 eV.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas
  Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import numpy as np

hbar = 1.0545718e-34      # Planck's constant over 2pi (J.s)
e = 1.60217662e-19        # elementary charge (C)
m_e = 9.10938356e-31      # electron mass (kg)
epsilon_0 = 8.854187817e-12  # vacuum permittivity (F/m)
c = 3e8                   # speed of light (m/s)

# ~~~~~~~~~~~~~~~~: Inputs :~~~~~~~~~~~~~~~~
lambda_nm = 1064                   # wavelength in nm
I0 = 5e13 * 1e4                    # convert W/cm^2 to W/m^2
Ip_eV = 13.6                       # Ionization potential in eV
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lambda_m = lambda_nm * 1e-9        # convert to meters
Ip_J = Ip_eV * e                   # in Joules

# Calculate laser frequency omega
omega = 2 * np.pi * c / lambda_m   # angular frequency (rad/s)

# Peak electric field amplitude E0
E0 = np.sqrt(2 * I0 / (c * epsilon_0))  # V/m

# Keldysh parameter
gamma = omega * np.sqrt(2 * m_e * Ip_J) / (e * E0)

print('gamma = ', gamma)

