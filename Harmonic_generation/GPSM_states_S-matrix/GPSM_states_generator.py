"""
File: GPSM_state_generator.py
Project: HHG-SaDAS
Code Description:
    - *** [MAIN GPSM STATES GENERATING CODE] ***
    - it automatically creates data directory if not previously created.
    - it checks, if for a particular file is already existing or not.
    - it runs and computes only if the calculation is unique.
    - The format of a generated file is shown below.

state_file = 'He_States_SAE-M1__l=0_nos=10_N=200_rmax=200_Lmap=20.xlsx'

+-----+-------------+----------+----------+----------+-----+----------+
| Row | r(x) (a.u.) |   A(1s)  |   A(2s)  |   A(3s)  | ... |  A(10s)  |
+-----+-------------+----------+----------+----------+-----+----------+
|  1  | 0.00166     | -0.00041 |  4.6E-05 | 0.000118 | ... | 8.64E-06 |
|  2  | 0.005566    |  0.001819| -0.00021 | -0.00053 | ... | -3.9E-05 |
|  3  | 0.011707    | -0.00454 |  0.000513| 0.001318 | ... | 9.63E-05 |
|  4  | 0.020085    |  0.00877 | -0.00099 | -0.00254 | ... | -0.00019 |
| ... | ...         | ...      | ...      | ...      | ... | ...      |
+-----+-------------+----------+----------+----------+-----+----------+
| 199 | 199.7993166 |-7.14E-13 |-1.38E-10 | 5.80E-11 | ... |-0.000331 |
+-----+-------------+----------+----------+----------+-----+----------+


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
total_states = 10
conf_info_string = conf_pot_selector(confinement_model, 0)[1]

if not confined:
    file_name = f'{evolving_atom}_States_{SAE_model}__l={l}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.xlsx'
else: file_name = f'{evolving_atom}@C60_States_{SAE_model}__l={l}_{conf_info_string}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.xlsx'

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
    data_wavefunction[f'A({state_name(Eth + l + 1, l)})'] = A[Eth]
df_A_r = pd.DataFrame(data_wavefunction)
df_A_r.to_excel(file_path, index=False)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Printing eigenvalues :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(f'Total number of states  :', total_states)
[print(f'E[{i}]~ {state_name(i + l + 1, l)} : {np.round(E[i], 10):.9f} a.u (Hartree)') for i in range(total_states)]

print(f"\npsi_file = '{file_name}'\n")
end_time = time.time()
print(f'Execution Time           : {end_time - start_time:.2f} seconds')
