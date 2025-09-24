"""
File: time_evolution.py
Project: HHG-SaDAS
Code Description:
    | *** [Main time evolution code] ***
    | Evolves initial m >= 0 states for n time steps using partial-wave expansion.
    | Computes electric field and dipole moment, saving the data in Excel.
    | Wavefunction evolution is calculated only at radial and angular collocation points (ri, theta_j).

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Only partial waves are evolved: glm(t+dt/2) = S(l) * glm(t)
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))   # Ensure project root (HHG-SaDAS) is in sys.path

import time
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from Assistant.Time_conversion import secs_to_hr_min_sec
from Assistant.Decorate_axes import decorate_axes_L as da
from Harmonic_generation.parameters_and_functions import (
    n, l, m,                                                                                     # initial state
    t, roots, colloc_pt, theta_k,                                                                         # arrays
    N, L, r_max, L_map, k_max, l_max, r0, dt, time_step, evolving_atom, eta_t,                   # parameters
    f, g_lm, Y_lm, Absorber_func, state_name, E_field, V_int, dipole_moment, Ps, show_E_field    # functions
)


# ~~~~~~~~~~~~~~~~~~~~~~: Common Figure Settings :~~~~~~~~~~~~~~~~~~~~~
width = 6.2                         # Width in inches
height = 4                          # Height in inches
fig_scale_factor = 2                # big=2 ; medium=1.5; small=1
tickslabel_size = 18
label_fontsize = 19
fig_size = (fig_scale_factor*width, fig_scale_factor*height)

dec_color = np.concatenate((da.mc.C_L, da.mc.des_col_1))
plt.rc('font', **{'family': 'serif'})
plt.rcParams['axes.prop_cycle'] = da.cycler(color=dec_color)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~: Importing files :~~~~~~~~~~~~~~~~~~~~~~~
this_dir = Path(__file__).resolve().parent

"""
Specify the data file names for the 'compatible' GPSM_states and S_matrix.

[Extremely Important Notice]:
By "compatible," the parameters used in these files must match those defined 
in parameters_and_functions.py.

Example:
If parameters_and_functions.py defines:
    L_map = 80
but the chosen files are :
    psi_file = 'He_States_SAE-M1__l=0_nos=10_N=200_rmax=200_Lmap=20.xlsx'
    S_matrix_file = 'He_Smatrix_SAE-M1__m=0_lmax=20_kmax=50_N=200_r_max=200_L_map=20_dt=0.1.xlsx'
then the code will give wrong results. This is because the imported nonlinear radial mapping 
function in this script from parameters_and_functions.py:

    def f(x, Lmap=L_map):
        r"
        Nonlinear radial mapping function.
        ...
        "

produces a radial grid that does not match the one encoded in the data files.

In short: ensure that the parameters in the data file names are consistent 
with those in parameters_and_functions.py, otherwise the grid and data 
will be incompatible.

[For crosschecking]: running this file will show what parameter values 
parameters_and_functions.py (and this script) is currently using.
Make sure these parameters are matching with the given datafile names: `psi_file' and `S_matrix_file'.
"""

psi_file = 'He_States_SAE-M1__l=1_nos=10_N=200_rmax=200_Lmap=20.xlsx'
psi_file = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Free_atom' / psi_file
psi_data = pd.read_excel(psi_file, header=None, skiprows=1).to_numpy().T

S_matrix_file = 'He_Smatrix_SAE-M1__m=1_lmax=20_kmax=50_N=200_r_max=200_L_map=20_dt=0.1.xlsx'
S_matrix_full_path = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Free_atom' / S_matrix_file
S_matrix_data = pd.read_excel(S_matrix_full_path, header=None, skiprows=1).to_numpy().T
S_matrix = np.array([[complex(*map(float, elem.split(','))) for elem in column] for column in S_matrix_data]).reshape(l_max+1, N-1, N-1)



# ~~~~~~~: Some pre-computed arrays to make calculations faster :~~~~~~~
r = f(colloc_pt)                                      # Radial coordinate in a.u
A_r = psi_data[1:][n-1]                               # Being the eigenstate of matrix hamiltonian, we'll evolve A(r).
absorber = np.array([Absorber_func(ri) for ri in r])

len_r = len(r); len_k = len(theta_k)
Y_lm_cos_theta_j = np.array([[Y_lm(l_ind+m, m, roots[j]) for j in range(len_k)] for l_ind in range(l_max+1)])


init_gl = np.zeros((l_max+1, len_r), dtype=np.complex128)
init_gl[l-m] = A_r.astype(np.complex128)


print('~~~~~~~~~~~~~: Time Evolution :~~~~~~~~~~~~~')
print('Evolving atom            :', evolving_atom)
print(f'Evolving initial state   : (n, l, m) : ({n+l}, {l}, {m}) ~', state_name(n + l, l))
print(f'θ_k[0]                   : {np.round(theta_k[0] * 180/np.pi, 4)} deg')
print(f'θ_k[-1]                  : {np.round(theta_k[-1] * 180/np.pi, 4)} deg')
print('S_matrix file name       :', S_matrix_file)
print('Initial state file name  :', psi_file)
print('Absorber radius (r_0)    :', r0)
print('Total time steps         :', time_step)
print('Estimated time (h, m, s) :', secs_to_hr_min_sec(eta_t * time_step), '\n')

if show_E_field:
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    da.decorate_2d(ax1)

    E_array = [E_field(ti) for ti in t]
    ax1.plot(t, E_array)
    ax1.set_title(f'E(t) ~ max step = {time_step}', fontsize=15)
    ax1.axvline(t[time_step], color='red', label=f't[{time_step}]')
    ax1.set_ylim(2*min(E_array), 2*max(E_array))
    ax1.fill_between(t, E_array, alpha=0.2)
    ax1.legend(loc='upper right', fontsize=15, framealpha=0.5, edgecolor='w')
    plt.show(block=False)       # Show plot without blocking code execution
    plt.pause(5)                # Pause for 5 seconds
    plt.close()                 # Close the plot

d_t_array = np.array([])          # Dipole moment array: d(t). Doesn't include initial wavefunction's dipole moment
population_den_array = np.array([])

dipole_moment_data = {'t (a.u)' : t[0:time_step],
                      'E(t)'    : [E_field(ti) for ti in t[0:time_step]]}
zero_psi = np.zeros((len_k, len_r), dtype=np.complex128); gl_psi = zero_psi
start_time = time.time()

gl_empty = np.empty((l_max+1, len(r)), dtype=np.complex128)         # Empty gl_array to be passed in gl() function.
for ti in range(time_step):
    psi_1 = 0 * zero_psi             # ψ1(r, θ) = exp{-iH0(dt)/2} • ψ0(r, θ)
    psi_2 = 0 * zero_psi             # ψ2(r, θ) = exp{-iV(r, θ, t+dt/2)(dt)/2} • ψ1(r, θ)

    for j in range(len_k):
        for l_index in range(l_max):
            psi_1[j] += np.dot(S_matrix[l_index], init_gl[l_index]) * Y_lm_cos_theta_j[l_index, j]
        psi_2[j] = np.exp(-1j * V_int(r, theta_k[j], t[ti]+dt/2) * dt) * psi_1[j]

    gl_2_array = g_lm(psi_2, gl_empty)
    for l_index in range(l_max):
        init_gl[l_index] = np.dot(S_matrix[l_index], gl_2_array[l_index]) * absorber

    print(f'Evolution step {ti:<5}: {((ti + 1) / time_step) * 100:.4f}%')

    dipole_mom = dipole_moment(r, init_gl)
    population_den_array = np.append(population_den_array, Ps(init_gl))
    d_t_array = np.append(d_t_array, dipole_mom)

end_time = time.time()


# ~~~~~~~~~~~~~~~~~~~: Saving dipole moment :~~~~~~~~~~~~~~~~~~~
# In the dipole moment files where cpp is not mentioned, assumed to be cpp=60
dipole_file_name = f'd_len_Ps_{time_step}_{evolving_atom}_{state_name(n+l, l)}_m={m}__L={L}_k_max={k_max}_N={N}_r_max={r_max}_L_map={L_map}_dt={dt}_wAb_r0={r0}_all_ok.xlsx'
dipole_moment_data['d(t)'] = d_t_array
dipole_moment_data['Ps(t)'] = population_den_array
df_dipole_moment_data = pd.DataFrame(dipole_moment_data)
d_file_path = rf'E:\Python_programs\HHG\GPSM\GPSM_Y_lm\free_SAE\Important_files\{dipole_file_name}'
# df_dipole_moment_data.to_excel(d_file_path, index=False)
print(f"'{dipole_file_name}'")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: plotting :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig2 = plt.figure(figsize=(18, 9))
ax2 = fig2.add_subplot(311)                         # Electric Field
ax3 = fig2.add_subplot(312)                         # dipole moment, from t=0
ax4 = fig2.add_subplot(313)                         # dipole moment, from t=0
da.decorate_2d([ax2, ax3, ax4])

ax3.plot(d_t_array, lw=2, color='deeppink', label='d(t)')
ax2.plot([E_field(ti) for ti in t[0:time_step]], lw=2, color='#58C4DD', label='E(t)')
ax4.plot(population_den_array, lw=2, color='orangered', label=r'P$_s$(t)')

ax3.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
ax2.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
ax4.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
fig2.suptitle(dipole_file_name, fontsize=13)

fig2.subplots_adjust(
    top=0.92,
    bottom=0.06,
    left=0.048,
    right=0.97,
    hspace=0.275,
    wspace=0.205
)

print('\n')
print('Average time for each step     :', (end_time-start_time)/time_step, ' sec')
print('Total Execution Time (h, m, s)  :', secs_to_hr_min_sec(end_time - start_time))

plt.show()
