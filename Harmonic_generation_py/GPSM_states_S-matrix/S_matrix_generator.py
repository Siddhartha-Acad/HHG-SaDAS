"""
File: S_matrix_generator.py
Project: HHG-SaDAS
Code Description:
    - *** [MAIN S-MATRIX GENERATING CODE] ***
    - Calculates the S-matrix for different angular momentum states (l).
    - Summation over all states indexed by k_max (number of eigenstates S(l)-matrix is made of).
    - Outputs the S-matrix data for all 'l' values in a numpy binary format '.npy'


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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import time
from scipy.linalg import eigh
from Assistant.Time_conversion import secs_to_hr_min_sec
from Harmonic_generation_py.parameters_and_functions import *


# ~~~~~~~~~~~~~~~~~~~~~~~~: File name and data arrangement system :~~~~~~~~~~~~~~~~~~~~~~~~
conf_info_string = conf_selector(confinement_model, 0)[1]

if not confined:
    file_name = f'{evolving_atom}_Smatrix_{SAE_model}_m={m}_lmax={l_max}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.npy'
else: file_name = f'{evolving_atom}@C60_Smatrix_{SAE_model}_m={m}_{conf_info_string}_lmax={l_max}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.npy'

if confined:
    output_dir = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_S-matrix_data' / 'Confined_atom'
else:
    output_dir = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_S-matrix_data' / 'Free_atom'
output_dir.mkdir(parents=True, exist_ok=True)  # Create if it doesn't exist

file_path = output_dir / file_name
if file_path.exists():
    raise FileExistsError(f"File already exists in dir: '{output_dir}'\n'{file_name}'")


print('Azimuthal quantum num. (m)  :', m)
print(f'S matrix range              : S({m}) to S({m+l_max})')
print('total S matrix (l_max+1)    :', l_max+1, '\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Computing S-matrix :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for l in range(m, l_max+m+1):
    # ~~~~~~~~~~~~~~~~~~~~~: H matrix, Eigenvalues (E), Eigenvectors (A) :~~~~~~~~~~~~~~~~~~~~~
    H_matrix = np.zeros((N - 1, N - 1))

    for i in range(N - 1):
        for j in range(i, N - 1):                   # Only computing the upper triangle
            H_matrix[i, j] = H(l, i, j, model=SAE_model)

    H_matrix += np.triu(H_matrix, 1).T              # Mirror to lower triangle (excluding diagonal)

    E, A = eigh(H_matrix, subset_by_index=[0, k_max-1])
    A = A.T


    S_matrix = (A.T * np.exp(-1j * E * dt / 2)) @ A                     # A: (k_max, n); A.T*(phase): (n, k_max)
    
    print(l, np.real(S_matrix[0, 0:5]))

