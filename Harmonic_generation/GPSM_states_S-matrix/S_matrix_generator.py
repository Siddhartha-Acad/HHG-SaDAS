"""
File: S_matrix_generator.py
Project: HHG-SaDAS
Code Description:
    - *** [MAIN S-MATRIX GENERATING CODE] ***
    - Calculates the S-matrix for different angular momentum states (l).
    - Summation over all states indexed by k_max (number of eigenstates S(l)-matrix is made of).
    - Outputs a file where each column corresponds to a specific l value,
      containing the flattened complex S_matrix.
    - Each matrix element is stored as a string in the format: "Real_part, Imaginary_part".

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
import pandas as pd
from scipy.linalg import eigh
from Harmonic_generation.parameters_and_functions import *
start_time = time.time()


# ~~~~~~~~~~~~~~~~~~~~~~~~: File name and data arrangement system :~~~~~~~~~~~~~~~~~~~~~~~~
conf_info_string = conf_pot_selector(confinement_model, 0)[1]

if not confined:
    file_name = f'{evolving_atom}_Smatrix_{SAE_model}__m={m}_lmax={l_max}_kmax={k_max}_N={N}_r_max={r_max}_L_map={L_map}_dt={dt}.xlsx'
else: file_name = f'{evolving_atom}@C60_Smatrix_{SAE_model}__m={m}_{conf_info_string}_lmax={l_max}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.xlsx'

if confined:
    output_dir = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Confined_atom'
else:
    output_dir = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Free_atom'
output_dir.mkdir(parents=True, exist_ok=True)  # Create if it doesn't exist

file_path = output_dir / file_name
if file_path.exists():
    raise FileExistsError(f"File already exists in dir: '{output_dir}'\n'{file_name}'")


print('Azimuthal quantum num. (m)  :', m)
print(f'S matrix range              : S({m}) to S({m+l_max})')
print('total S matrix (l_max+1)    :', l_max+1, '\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Computing S-matrix :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_S_matrix = {}
energy_eigenvalues = {}
for l in range(m, l_max+m+1):
    # ~~~~~~~~~~~~~~~~~~~~~: H matrix, Eigenvalues (E), Eigenvectors (A) :~~~~~~~~~~~~~~~~~~~~~
    H_matrix = np.zeros((N - 1, N - 1))
    for i in range(N - 1):
        for j in range(i, N - 1):                   # Only computing the upper triangle
            H_matrix[i, j] = H_matrix[j, i] = H(l, i, j, model=SAE_model)

    E, A = eigh(H_matrix, subset_by_index=[0, k_max-1])
    A = A.T

    energy_eigenvalues[f'l={l}'] = E                # Store eigenvalues for l
    print(f'S-matrix for l={l:<2}:  DONE')
    # positive_energy_states = np.sum(E > 0)
    # negative_energy_states = np.sum(E < 0)
    # print(f'negative energy states (E<0) : {negative_energy_states}')
    # print(f'positive energy states (E>0) : {positive_energy_states}\n')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: S matrix :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    S_matrix_real = np.zeros((N-1, N-1))
    S_matrix_imag = np.zeros((N-1, N-1))
    for i in range(N - 1):
        for j in range(i, N - 1):                               # Only compute the upper triangle
            matrix_ele = S(E, A, i, j)
            S_matrix_real[i][j] = np.real(matrix_ele)
            S_matrix_imag[i][j] = np.imag(matrix_ele)
            if i != j:                                          # Mirror the values to the lower triangle
                S_matrix_real[j][i] = np.real(matrix_ele)
                S_matrix_imag[j][i] = np.imag(matrix_ele)
    flat_S_matrix_real = S_matrix_real.flatten()
    flat_S_matrix_imag = S_matrix_imag.flatten()
    data_S_matrix[f'l={l}'] = [f'{flat_S_matrix_real[i]}, {flat_S_matrix_imag[i]}' for i in range((N - 1) * (N - 1))]


# ~~~~~~~~~~~~~~~~~~~~~~~~: Writing S-matrices data to .xlsx file :~~~~~~~~~~~~~~~~~~~~~~~~
df_S_matrix = pd.DataFrame(data_S_matrix)
df_S_matrix.to_excel(file_path, index=False)
print(f"\nS_matrix_file = '{file_name}'\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: saving EgVals: .txt :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if save_Egvals_with_Smatrix:
    if not confined:
        output_name = f'{evolving_atom}_EgVals__m={m}_lmax={l_max}_N={N}_rmax={r_max}_Lmap={L_map}.txt'
    else:
        output_name = f'{evolving_atom}@C60_EgVals__m={m}_{conf_info_string}_lmax={l_max}_N={N}_rmax={r_max}_Lmap={L_map}.txt'
    output_path = output_dir / output_name

    with open(output_path, 'w') as f:
        f.write(" ".join([f"l={l}" for l in range(m, l_max+m+1)]) + "\n")
        max_rows = max(len(vals) for vals in energy_eigenvalues.values())
        for row in range(max_rows):
            row_data = []
            for l in range(m, l_max+m+1):
                row_data.append(f"{energy_eigenvalues[f'l={l}'][row]:.6f}" if row < len(energy_eigenvalues[f'l={l}']) else "")
            f.write(" ".join(row_data) + "\n")
    print(f"EgVals saved : '{output_name}'")

end_time = time.time()
print(f'Execution Time : {end_time - start_time:.2f} seconds')
