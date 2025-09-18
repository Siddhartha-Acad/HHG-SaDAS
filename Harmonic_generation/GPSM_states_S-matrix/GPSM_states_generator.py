"""
~ GPSM_state_generator.py

• r = f(x): nonlinear mapped radial coordinate.
• Takes l value from H_ij_gen_anim.py and calculates eigenstates A(r).
• Calculates φ(r), psi(r) = u(r).
• Generates file containing all A(r) corresponding to l.
• These file is used to get the initial wavefunction in Evolve_dt.py, n_evolution.py, save_frame.py.
"""

import time
import pandas as pd
from scipy.linalg import eigh
from Harmonic_generation.parameters_and_functions import *
start_time = time.time()


conf_info_string = conf_pot_selector(confinement_model, 0)[1]
# ~~~~~~~~~~~~~~~~~~~~~: H matrix, Eigenvalues (E), Eigenvectors (A) :~~~~~~~~~~~~~~~~~~~~~
H_matrix = np.zeros((N - 1, N - 1))
for i in range(N - 1):
    for j in range(i, N - 1):                           # Only computing the upper triangle
        H_matrix[i, j] = H_matrix[j, i] = H(l, i, j, model=SAE_model)

total_states = 10
E, A = eigh(H_matrix, subset_by_index=[0, total_states-1])
A = A.T


if not confined:
    file_name = f'{evolving_atom}_States__l={l}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.xlsx'
else: file_name = f'{evolving_atom}@C60_States__l={l}_{conf_info_string}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.xlsx'

# ~~~~~~~~~~~~~~~~~~~~~~~~: Writing GPSM-states data to .xlsx file :~~~~~~~~~~~~~~~~~~~~~~~
r = f(colloc_pt)                                    # radial coordinate in atomic unit. (Nonlinearly discretised)
data_wavefunction = {'r (a.u.)': r}                 # First column of the data file is the radial grid.
for Eth in range(total_states):
    data_wavefunction[f'A_{Eth}'] = A[Eth]
df_A_r = pd.DataFrame(data_wavefunction)
df_A_r.to_excel(file_name, index=False)


print(f'Total number of states  :', total_states)
[print(f'E[{i}]~ {state_name(i + 1 + l, l)} : {np.round(E[i], 10):.9f} a.u (Hartree)') for i in range(total_states)]

print(f"\n'{file_name}'\n")
end_time = time.time()
print(f'Execution Time           : {end_time - start_time:.2f} seconds')
