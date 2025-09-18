"""
File: GPSM_state_generator.py
Project: HHG-SaDAS
Code Description:
    - *** [MAIN GPSM STATES GENERATING CODE] ***
    - it automatically creates data directory if not previously created.
    - it checks, if for a particular file is already existing or not.
    - it runs and computes only if the calculation is unique.

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

import time
import pandas as pd
from scipy.linalg import eigh
from Harmonic_generation.parameters_and_functions import *
start_time = time.time()


# ~~~~~~~~~~~~~~~~~~~~~~~~: File name and data arrangement system :~~~~~~~~~~~~~~~~~~~~~~~~
total_states = 10
conf_info_string = conf_pot_selector(confinement_model, 0)[1]

if not confined:
    file_name = f'{evolving_atom}_States__l={l}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.xlsx'
else: file_name = f'{evolving_atom}@C60_States__l={l}_{conf_info_string}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.xlsx'

if confined:
    output_dir = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Confined_atom'
else:
    output_dir = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Free_atom'
output_dir.mkdir(parents=True, exist_ok=True)  # Create if it doesn't exist

file_path = output_dir / file_name
if file_path.exists():
    raise FileExistsError(f"File already exists:\nbound_states_file = '{file_name}'")


# ~~~~~~~~~~~~~~~~~~~~~: H matrix, Eigenvalues (E), Eigenvectors (A) :~~~~~~~~~~~~~~~~~~~~~
H_matrix = np.zeros((N - 1, N - 1))
for i in range(N - 1):
    for j in range(i, N - 1):                       # Only computing the upper triangle
        H_matrix[i, j] = H_matrix[j, i] = H(l, i, j, model=SAE_model)
E, A = eigh(H_matrix, subset_by_index=[0, total_states-1])
A = A.T


# ~~~~~~~~~~~~~~~~~~~~~~~~: Writing GPSM-states data to .xlsx file :~~~~~~~~~~~~~~~~~~~~~~~
r = f(colloc_pt)                                    # radial coordinate in atomic unit. (Nonlinearly discretised)
data_wavefunction = {'r (a.u.)': r}                 # First column of the data file is the radial grid.
for Eth in range(total_states):
    data_wavefunction[f'A_{Eth}'] = A[Eth]
df_A_r = pd.DataFrame(data_wavefunction)
df_A_r.to_excel(file_path, index=False)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Printing eigenvalues :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(f'Total number of states  :', total_states)
[print(f'E[{i}]~ {state_name(i + 1 + l, l)} : {np.round(E[i], 10):.9f} a.u (Hartree)') for i in range(total_states)]

print(f"\n'{file_name}'\n")
end_time = time.time()
print(f'Execution Time           : {end_time - start_time:.2f} seconds')
