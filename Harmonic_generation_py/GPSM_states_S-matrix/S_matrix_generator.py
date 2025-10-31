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
import argparse
from scipy.linalg import eigh
from Assistant.Time_conversion import secs_to_hr_min_sec
from Harmonic_generation_py.parameters import *
from Harmonic_generation_py.functions import *

parser = argparse.ArgumentParser()
parser.add_argument("-v", action="store_true")
args = parser.parse_args()

if args.v:
    print_info()

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
    print(f"File already exists : {file_path.name}\n")
    sys.exit(0)                         # Exit program gracefully


print('Azimuthal quantum num. (m)  :', m)
print(f'S matrix range              : S({m}) to S({m+l_max})')
print('total S matrix (l_max+1)    :', l_max+1, '\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Computing S-matrix :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_S_matrix = []
energy_eigenvalues = {}
start_time = time.perf_counter()

for l in range(m, l_max+m+1):
    # ~~~~~~~~~~~~~~~~~~~~~: H matrix, Eigenvalues (E), Eigenvectors (A) :~~~~~~~~~~~~~~~~~~~~~
    H_matrix = np.zeros((N - 1, N - 1))

    for i in range(N - 1):
        for j in range(i, N - 1):                   # Only computing the upper triangle
            H_matrix[i, j] = H(l, i, j, model=SAE_model)

    H_matrix += np.triu(H_matrix, 1).T              # Mirror to lower triangle (excluding diagonal)

    E, A = eigh(H_matrix, subset_by_index=[0, k_max-1])
    A = A.T

    energy_eigenvalues[f'l={l}'] = E                # Store eigenvalues for l
    print(f'S-matrix for l={l:<3}:  DONE')
    # positive_energy_states = np.sum(E > 0)
    # negative_energy_states = np.sum(E < 0)
    # print(f'negative energy states (E<0) : {negative_energy_states}')
    # print(f'positive energy states (E>0) : {positive_energy_states}\n')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: S matrix :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    S_matrix = (A.T * np.exp(-1j * E * dt / 2)) @ A                     # A: (k_max, n); A.T*(phase): (n, k_max)
    
    data_S_matrix.append(S_matrix)


# ~~~~~~~~~~~~~~~~~~~~~~~~: Writing S-matrices data to .npy file :~~~~~~~~~~~~~~~~~~~~~~~~
data_S_matrix = np.array(data_S_matrix, dtype=np.complex128)            # shape: (l_max+1, N-1, N-1)
np.save(file_path, data_S_matrix)
print(f"\nS_matrix_file = '{file_name}'")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: saving EgVals: .txt :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if save_Egvals_with_Smatrix:
    if not confined:
        output_name = f'{evolving_atom}_EgVals_lmax={l_max}_N={N}_rmax={r_max}_Lmap={L_map}.txt'
    else:
        output_name = f'{evolving_atom}@C60_EgVals_{conf_info_string}_lmax={l_max}_N={N}_rmax={r_max}_Lmap={L_map}.txt'

    output_path = output_dir / output_name

    if output_path.exists():
        print(f"File already exists: '{output_name}' — skipping.")
    else:
        with open(output_path, 'w') as f:
            f.write(" ".join([f"l={l}" for l in range(m, l_max+m+1)]) + "\n")
            max_rows = max(len(vals) for vals in energy_eigenvalues.values())
            for row in range(max_rows):
                row_data = []
                for l in range(m, l_max+m+1):
                    row_data.append(
                        f"{energy_eigenvalues[f'l={l}'][row]:.6f}"
                        if row < len(energy_eigenvalues[f'l={l}']) else ""
                    )
                f.write(" ".join(row_data) + "\n")
        print(f"EgVals_file = '{output_name}'")


end_time = time.perf_counter()
wall_time = end_time - start_time

if wall_time > 300.0:
    print(f'\nExecution Wall-Time (h, m, s) : {secs_to_hr_min_sec(wall_time)}')
else:
    print(f'\nExecution Wall-Time : {wall_time:.3f} seconds')
