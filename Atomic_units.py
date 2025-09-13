"""
File: Atomic_units.py
Project: HHG-SaDAS
Code Description:
    | Utility module for defining and handling physical constants
    | (SI and atomic units) used in simulations of High Harmonic
    | Generation (HHG) from free and confined atoms.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- All constants follow CODATA 2018 recommended values unless stated otherwise.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

from scipy.constants import c, epsilon_0, mu_0, physical_constants


def omega_au(lambda_nm):
    """
    lambda_nm : wavelength in nm.
    returns : angular frequency in atomic unit.
    """
    return 45.5633526 / lambda_nm

def wavelength_nm(w_au):
    """
    omega_au : angular frequency in atomic unit.
    returns : wavelength in nm.
    """
    return 45.5633526 / w_au


def I_to_E_0(I_w_per_cm_sq):
    I = I_w_per_cm_sq * 10**4               # converting to W/m^2
    E_0 = np.sqrt(2 * I / (c * epsilon_0))
    E_0_au = E_0 / E_fld_0
    return E_0_au



# Fundamental constants
c         = c                                                        # Speed of light (m/s)
epsilon_0 = epsilon_0                                                # Free space permittivity (F/m)

# Atomic units (CODATA values)
a0        = physical_constants['Bohr radius'][0]                     # Length (m)    : Bohr radius
v0        = physical_constants['atomic unit of velocity'][0]         # Velocity (m/s): Velocity of an electron in the first Bohr orbit
T0        = physical_constants['atomic unit of time'][0]             # Time (s)      : Time for an electron to travel distance a0 at velocity v0
nu0       = 1 / T0                                                   # Frequency (Hz): Inverse of the unit of time
Energy_0  = physical_constants['Hartree energy in eV'][0]            # Energy (eV)   : Twice the binding energy of hydrogen
E_fld_0   = physical_constants['atomic unit of electric field'][0]   # Electric field strength (V/m)
Int_0     = (E_fld_0**2) / (2 * mu_0 * c) * 1e-4                     # Energy flux (intensity) (W/cm²)


if __name__ == '__main__':
    import numpy as np

    print(I_to_E_0(5 * 10**13))
    print(np.sqrt(5*10**13 / Int_0))