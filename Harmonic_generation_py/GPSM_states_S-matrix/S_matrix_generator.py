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
    - An example of a generated S-matrix file is shown below.

S_matrix_file = 'He_Smatrix_SAE-M1__m=1_lmax=20_kmax=50_N=200_r_max=200_L_map=80_dt=0.1.xlsx'

+-----+---------------+---------------+---------------+---------------+-----+---------------+
| Row |    S(l=1)     |    S(l=2)     |    S(l=3)     |    S(l=4)     | ... |   S(l=21)     |
+-----+---------------+---------------+---------------+---------------+-----+---------------+
|  1  | (1.028e-14,   | (4.204e-22,   | (3.431e-24,   | (2.106e-24,   | ... | (2.902e-24,   |
|     |  -3.338e-17)  |  -3.592e-24)  |  -1.813e-26)  |  -1.942e-26)  |     |  -5.925e-26)  |
+-----+---------------+---------------+---------------+---------------+-----+---------------+
|  2  | (-1.545e-13,  | (-2.159e-20,  | (2.993e-23,   | (1.965e-23,   | ... | (3.212e-23,   |
|     |  5.016e-16)   |  1.953e-22)   |  -1.612e-25)  |  -1.812e-25)  |     |  -6.558e-25)  |
+-----+---------------+---------------+---------------+---------------+-----+---------------+
|  3  | (8.164e-13,   | (2.398e-19,   | (9.406e-23,   | (6.743e-23,   | ... | (1.387e-22,   |
|     |  -2.651e-15)  |  -2.149e-21)  |  -4.245e-25)  |  -6.211e-25)  |     |  -2.833e-24)  |
+-----+---------------+---------------+---------------+---------------+-----+---------------+
|  4  | (-2.726e-12,  | (-1.379e-18,  | (2.017e-22,   | (1.492e-22,   | ... | (3.957e-22,   |
|     |  8.852e-15)   |  1.239e-20)   |  -1.771e-24)  |  -1.374e-24)  |     |  -8.080e-24)  |
+-----+---------------+---------------+---------------+---------------+-----+---------------+
|  5  | (7.009e-12,   | (5.436e-18,   | (2.933e-22,   | (2.633e-22,   | ... | (8.895e-22,   |
|     |  -2.276e-14)  |  -4.879e-20)  |  2.741e-24)   |  -2.413e-24)  |     |  -1.816e-23)  |
+-----+---------------+---------------+---------------+---------------+-----+---------------+
|  6  | (-1.521e-11,  | (-1.681e-17,  | (6.506e-22,   | (4.020e-22,   | ... | (1.712e-21,   |
|     |  4.939e-14)   |  1.510e-19)   |  -2.223e-23)  |  -3.737e-24)  |     |  -3.496e-23)  |
+-----+---------------+---------------+---------------+---------------+-----+---------------+
| ... | ...           | ...           | ...           | ...           | ... | ...           |
+-----+---------------+---------------+---------------+---------------+-----+---------------+
|39601| (1.311e-03,   | (1.360e-03,   | (1.422e-03,   | (1.464e-03,   | ... | (2.167e-03,   |
|     |  -1.134e-05)  |  -1.207e-05)  |  -1.302e-05)  |  -1.368e-05)  |     |  -2.703e-05)  |
+-----+---------------+---------------+---------------+---------------+-----+---------------+


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
    file_name = f'{evolving_atom}_Smatrix_{SAE_model}__m={m}_lmax={l_max}_kmax={k_max}_N={N}_r_max={r_max}_L_map={L_map}_dt={dt}.npy'
else: file_name = f'{evolving_atom}@C60_Smatrix_{SAE_model}__m={m}_{conf_info_string}_lmax={l_max}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.npy'

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Computing S-matrix :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_S_matrix = []
energy_eigenvalues = {}
start_time = time.perf_counter()

for l in range(m, l_max+m+1):
    # ~~~~~~~~~~~~~~~~~~~~~: H matrix, Eigenvalues (E), Eigenvectors (A) :~~~~~~~~~~~~~~~~~~~~~
    H_matrix = np.zeros((N - 1, N - 1))
    for i in range(N - 1):
        for j in range(i, N - 1):                   # Only computing the upper triangle
            H_matrix[i, j] = H_matrix[j, i] = H(l, i, j, model=SAE_model)

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


# ~~~~~~~~~~~~~~~~~~~~~~~~: Writing S-matrices data to .npz file :~~~~~~~~~~~~~~~~~~~~~~~~
data_S_matrix = np.array(data_S_matrix, dtype=np.complex128)            # shape: (l_max+1, N-1, N-1)
np.save(file_path, data_S_matrix)
print(f"\nS_matrix_file = '{file_name}'")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: saving EgVals: .txt :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if save_Egvals_with_Smatrix:
    if not confined:
        output_name = f'{evolving_atom}_EgVals__lmax={l_max}_N={N}_rmax={r_max}_Lmap={L_map}.txt'
    else:
        output_name = f'{evolving_atom}@C60_EgVals__{conf_info_string}_lmax={l_max}_N={N}_rmax={r_max}_Lmap={L_map}.txt'

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
    print(f'\nWall Time (h, m, s) : {secs_to_hr_min_sec(wall_time)}')
else:
    print(f'\nWall Time : {wall_time:.3f} seconds')
