"""
File: GPSM_state_generator.py
Project: HHG-SaDAS
Code Description:
    - *** [MAIN GPSM STATES GENERATING CODE] ***
    - it automatically creates data directory if not previously created.
    - it checks, if for a particular file is already existing or not.
    - it runs and computes only if the calculation is unique.
    - The format of a generated file is shown below.

state_file = 'He_States_SAE-M1__l=0_nos=10_N=200_rmax=200_Lmap=20.dat'

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
import argparse
from scipy.linalg import eigh
from Harmonic_generation_py.parameters import *
from Harmonic_generation_py.functions import print_info, conf_selector, state_name, colloc_pt, H, f
start_time = time.perf_counter()

parser = argparse.ArgumentParser()
parser.add_argument("-v", action="store_true")
args = parser.parse_args()

if args.v:
    print_info()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     File name and data arrangement system      |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if not confined:
    file_name = f'{evolving_atom}_States_{SAE_model}_l={l}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.dat'
else: file_name = f'{evolving_atom}@C60_States_{SAE_model}_{confinement_model}_l={l}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.dat'

data_dir = 'Confined_atom' if confined else 'Free_atom'
file_path = this_dir / 'GPSM_states_S-matrix' / 'data_GPSM_states_S-matrix' / data_dir / file_name

if file_path.exists():
    print(f"File already exists : {file_path.name}\n")
    sys.exit(0)                         # Exit program gracefully


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  H matrix, Eigenvalues (E), Eigenvectors (A)   |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
H_matrix = np.zeros((N - 1, N - 1))

for i in range(N - 1):
    for j in range(i, N - 1):                   # Only computing the upper triangle
        H_matrix[i, j] = H(l, i, j, model=SAE_model)

H_matrix += np.triu(H_matrix, 1).T              # Mirror to lower triangle (excluding diagonal)

E, A = eigh(H_matrix, subset_by_index=[0, total_states-1])
A = A.T


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     Writing GPSM-states data to .dat file      |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r = f(colloc_pt)                                       # radial coordinate in atomic unit. (Nonlinearly discretised)

header = "r(a.u.) " + " ".join([f"A({state_name(Eth + l + 1, l)})" for Eth in range(total_states)])
data = np.column_stack([r] + [A[Eth] for Eth in range(total_states)])
np.savetxt(file_path, data, header=header, comments='', fmt='%.16e')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#              Printing eigenvalues              |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(f'Total number of states  :', total_states)
[print(f'E[{i}]~ {state_name(i + l + 1, l)} : {np.round(E[i], 10):.9f} a.u (Hartree)') for i in range(total_states)]

print(f"\nstate_file = '{file_name}'\n")

end_time = time.perf_counter()
print(f'Wall Time : {end_time - start_time:.3f} seconds')

