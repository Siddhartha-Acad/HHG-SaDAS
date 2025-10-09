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
from parameters_and_functions import (
    n, l, m,                                                                                                    # initial state
    t, roots, colloc_pt, theta_k,                                                                               # arrays
    show_E_field, print_serial_prog, confined,                                                                  # booleans
    f, g_lm, Y_lm, conf_selector, Absorber_func, state_name, E_field, V_int, dipole_moment, Ps,                 # functions
    N, L, r_max, L_map, k_max, l_max, r0, dt, time_step, evolving_atom, eta_t, SAE_model, confinement_model     # parameters
)

# ~~~~~~~~~~~~~~~~~~~~~~: Common Figure Settings :~~~~~~~~~~~~~~~~~~~~~
width = 6.2                         # Width in inches
height = 4                          # Height in inches
fig_scale_factor = 2                # big=2 ; medium=1.5; small=1
tickslabel_size = 18
label_fontsize = 19
fig_size = (fig_scale_factor*width, fig_scale_factor*height)

plt.rc('font', **{'family': 'serif'})
plt.rcParams['axes.prop_cycle'] = da.cycler(color=da.mc.C_L)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Importing files                           |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
this_dir = Path(__file__).resolve().parent

if confined:
    data_dir = 'Confined_atom'
else:
    data_dir = 'Free_atom'

"""
Specify the data file names for the 'compatible' GPSM_states and S_matrix.

[Extremely Important Notice]:
By "compatible," the parameters used in these files must match those defined 
in parameters_and_functions.py.

Example:
If parameters_and_functions.py defines:
    L_map = 80
but the chosen files are :
    state_file = 'He_States_SAE-M1__l=0_nos=10_N=200_rmax=200_Lmap=20.xlsx'
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
Make sure these parameters are matching with the given datafile names: `state_file' and `S_matrix_file'.
"""

state_file = 'H_States_SAE-M1__l=0_nos=10_N=200_rmax=200_Lmap=80.xlsx'
state_file = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / data_dir / state_file
state_data = pd.read_excel(state_file, header=None, skiprows=1).to_numpy().T

S_matrix_file = 'H_Smatrix_SAE-M1__m=0_lmax=20_kmax=50_N=200_r_max=200_L_map=80_dt=0.1.xlsx'
S_matrix_full_path = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / data_dir / S_matrix_file
S_matrix_data = pd.read_excel(S_matrix_full_path, header=None, skiprows=1).to_numpy().T
S_matrix = np.array([[complex(*map(float, elem.split(','))) for elem in column] for column in S_matrix_data]).reshape(l_max+1, N-1, N-1)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Printing system info. and show Electric field.           |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if confined:
    conf_string = conf_selector(confinement_model, 0)[1]

    print("~~~~~~~~~~~~~: Conf info :~~~~~~~~~~~~~")
    parts = conf_string.split('_')                  # Split by underscores
    conf_model = parts[0]                           # The first part (before first '_') is the confinement model
    params = [p for p in parts[1:] if '=' in p]     # Remaining parts contain key=value pairs

    print(f"{'Conf. model':<10}: {conf_model}")
    for p in params:
        key, val = p.split('=')
        print(f"{key:<10}: {val}")

print('~~~~~~~~~~~: Time Evolution :~~~~~~~~~~')
print('Evolving atom            :', evolving_atom)
print(f'Evolving initial state   : (n, l, m) : ({n+l}, {l}, {m}) ~', state_name(n + l, l))
print(f'θ_k[0]                   : {np.round(theta_k[0] * 180/np.pi, 4)} deg')
print(f'θ_k[-1]                  : {np.round(theta_k[-1] * 180/np.pi, 4)} deg')
print('S_matrix file name       :', S_matrix_file)
print('Initial state file name  :', state_file)
print('Absorber radius (r_0)    :', r0)
print('Total time steps         :', time_step)
print('Estimated time (h, m, s) :', secs_to_hr_min_sec(eta_t * time_step), '\n')

if show_E_field:
    fig1 = plt.figure(figsize=fig_size)
    ax1 = fig1.add_subplot(111)
    da.decorate_2d(ax1)

    E_array = [E_field(ti) for ti in t]
    ax1.plot(t, E_array)
    ax1.set_title(f'E(t) (a.u.) : max step = {time_step}', fontsize=15)
    ax1.axvline(t[time_step], linestyle='--', color='royalblue', label=f't[{time_step}]')
    ax1.set_ylim(2 * min(E_array), 2 * max(E_array))
    ax1.fill_between(t, E_array, alpha=0.2)
    ax1.legend(loc='upper right', fontsize=15, framealpha=0.5, edgecolor='w')
    plt.show(block=False)       # Show plot without blocking code execution
    plt.pause(5)                # Pause for 5 seconds
    plt.close()                 # Close the plot



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   Arrays to hold evolution data                    |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r = f(colloc_pt)                                        # Nonlinear Radial collocation grid (a.u)
A_r = state_data[1:][n - 1]                             # Selecting the initial state from state_data array and storing in A_r
absorber = np.array([Absorber_func(ri) for ri in r])    # The absorbing layer.

init_glm = np.zeros((l_max + 1, N - 1), dtype=np.complex128)    # Empty zero array to hold all the time evolving partial waves.
init_glm[l - m] = A_r.astype(np.complex128)     # The initial state is set as the only living partial wave. (be careful with the indexing here)
#                                               # Instead of making full 3d initial state and calculating partial waves, I can directly
#         [Initializing: partial waves]         # use the A_nl(r) as the initial partial wave.
#                                               # Details in: Section~2.3.8 `Multiple-step time evolution', see Figure 2.23(a)
glm_empty = np.empty((l_max+1, N-1), dtype=np.complex128)  # Empty gl_array to be passed in gl() function.

d_t_array = np.array([])                                      # Dipole moment array: d(t). Doesn't include initial wavefunction's dipole moment
population_den_array = np.array([])                           # Array to store Population density
zero_psi = np.zeros((L+1, N-1), dtype=np.complex128)    # To initiate the wavefunction A(ri, θj) before each loop.
Evolution_data = {'t (a.u)' : t[0:time_step],                            # The first column of the dipole moment file is reserved for time.
                  'E(t)'    : [E_field(ti) for ti in t[0:time_step]]}    # And the second column is reserved for electric field.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    STARTING MAIN TIME EVOLUTION                    |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start_time = time.process_time()     # Start measuring execution time (serial time evolution)
last_percent = -1                    # keeps track of the last printed percent

Y_lm_cos_theta_j = np.array([[Y_lm(l_ind+m, m, roots[j]) for j in range(L+1)] for l_ind in range(l_max+1)])    # precomputed Spherical Harmonics
for ti in range(time_step):
    psi_1 = 0 * zero_psi             # ψ1(r, θ) = exp{-iH0(dt)/2} • ψ0(r, θ)                  # See Eq.~2.84
    psi_2 = 0 * zero_psi             # ψ2(r, θ) = exp{-iV(r, θ, t+dt/2)(dt)/2} • ψ1(r, θ)     # See Eq.~2.86

    for j in range(L+1):                    # angular grid index j
        for l_index in range(l_max):        # summation index on 'l' of Eq.~2.85.
            psi_1[j] += np.dot(S_matrix[l_index], init_glm[l_index]) * Y_lm_cos_theta_j[l_index, j]     # calculating ψ1(r, θ)
        psi_2[j] = np.exp(-1j * V_int(r, theta_k[j], t[ti]+dt/2) * dt) * psi_1[j]                       # calculating ψ2(r, θ)

    glm_tilde = g_lm(psi_2, glm_empty)                                                   # the \tilde{g}_{\ell}(r, t) of Eq.~2.88
    for l_index in range(l_max):                                                         # Again summation index on 'l' of Eq.~2.85.
        init_glm[l_index] = np.dot(S_matrix[l_index], glm_tilde[l_index]) * absorber     # Implementing Eq.~2.89 and updating init_glm with absorbing function
    # main time evolution algorithm ends here...

    if print_serial_prog:
        percent = int(((ti + 1) / time_step) * 100)
        if percent > last_percent:                                      # update progress only after one percent.
            print(f"Evolution step {ti:<6}: {percent}%")                # It will show how much the process is completed.
            last_percent = percent

    d_t_array = np.append(d_t_array, dipole_moment(r, init_glm))                         # calculating and storing the dipole moment of this instant.
    population_den_array = np.append(population_den_array, Ps(init_glm))                 # calculating and storing the population density

end_time = time.process_time()             # ending time measurement.
CPU_time = end_time - start_time           # Total CPU time for computing the total time evolution.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Saving dipole moment; survival probability & correlation functions |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if not confined:
    Evo_data_file = f'Evo_steps={time_step}_{evolving_atom}({state_name(n + l, l)})_m={m}_{SAE_model}__L={L}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.xlsx'
else:
    Evo_data_file = f'Evo_steps={time_step}_{evolving_atom}({state_name(n + l, l)})@C60_m={m}_{SAE_model}_{confinement_model}__L={L}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.xlsx'

Evolution_data['d(t)'] = d_t_array                      # In the dipole moment files where cpp is not mentioned, assumed to be cpp=60
Evolution_data['Ps(t)'] = population_den_array
df_Evolution_data = pd.DataFrame(Evolution_data)

output_dir = this_dir / 'Time_evolution_data'
output_dir.mkdir(parents=True, exist_ok=True)               # Create if it doesn't exist
d_file_path = output_dir / f'{Evo_data_file}'
df_Evolution_data.to_excel(d_file_path, index=False)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              plotting                              |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig2 = plt.figure(figsize=(16, 9))
ax2 = fig2.add_subplot(221)                         # Electric Field
ax3 = fig2.add_subplot(223)                         # dipole moment
ax4 = fig2.add_subplot(122)                         # survival probability
da.decorate_2d([ax2, ax3, ax4])

ax3.plot(d_t_array, lw=2, color='deeppink', label='d(t)')
ax2.plot([E_field(ti) for ti in t[0:time_step]], lw=2, color='#58C4DD', label='E(t)')
ax4.plot(population_den_array, lw=2, color='orangered', label=r'P$_s$(t)')

ax3.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
ax2.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
ax4.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
fig2.suptitle(Evo_data_file, fontsize=13)

fig2.subplots_adjust(
    top=0.92,
    bottom=0.06,
    left=0.048,
    right=0.97,
    hspace=0.275,
    wspace=0.205
)

print('\n')
print(f"evo_data_file_name = '{Evo_data_file}'")
print('\n')

if CPU_time > 300.0:
    print(f"Average CPU time per step (eta_t)      : {CPU_time / time_step:.3f} seconds")
    print(f'Total CPU time for all steps (h, m, s) : {secs_to_hr_min_sec(CPU_time)}')
else:
    print(f"Average CPU time per step (eta_t) : {CPU_time / time_step:.3f} seconds")
    print(f'Total CPU time for all steps      : {CPU_time:.3f} seconds')

plt.show()
