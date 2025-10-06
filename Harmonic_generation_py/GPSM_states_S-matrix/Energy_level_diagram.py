"""
File: Energy_level_diagram.py
Project: HHG-SaDAS
Code Description:
    This code shows the energy level diagram for a particular atom (free/confined).
    The data file was generated as a byproduct while executing S_matrix_generator.py
    Whether this file will be created or not, is decided by the boolean parameter 'save_Egvals_with_Smatrix' in parameters_and_functions.pyv


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

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))   # Ensure project root (HHG-SaDAS) is in sys.path

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import Assistant.Decorate_axes.decorate_axes_L as da
this_dir = Path(__file__).resolve().parent

# ~~~~~~~~~~~~~~~~~~~~~~: Common Figure Settings :~~~~~~~~~~~~~~~~~~~~~
width = 4                           # Width in inches
height = 5                          # Height in inches
fig_scale_factor = 2                # big=2 ; medium=1.5; small=1
fig_size = (fig_scale_factor*width, fig_scale_factor*height)

plt.rc('font', **{'family': 'serif'})
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig = plt.figure(figsize=fig_size)
ax1 = fig.add_subplot(111)
da.decorate_2d(ax1)


file_name = 'He_EgVals__lmax=20_N=200_rmax=200_Lmap=80.txt'
file = this_dir / 'GPSM_states_and_Smatrix_data' / 'Free_atom' / f'{file_name}'
data = np.loadtxt(file, skiprows=1).T

l_max = 5; E_max = 0
l_values = np.arange(0, l_max, 1)
for l in l_values:
    energy_levels = data[l][data[l] < E_max]
    x_vals = np.full_like(energy_levels, l)
    ax1.hlines(energy_levels, x_vals - 0.3, x_vals + 0.3, colors='b', lw=2)



ax1.set_xlabel(r'Azimuthal Quantum Number ($\ell$) $\longrightarrow$', fontsize=15)
ax1.set_ylabel(r'Energy Eigenvalues (a.u) $\longrightarrow$', fontsize=15, labelpad=5)
ax1.set_title(f'{file_name}', fontsize=16, pad=30)
fig.suptitle('Energy Level Diagram', fontsize=17)
ax1.set_xticks(ticks=l_values, labels=l_values)
fig.subplots_adjust(
    top=0.88,
    bottom=0.11,
    left=0.130,
    right=0.945,
    hspace=0.2,
    wspace=0.2
)
plt.show()
