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
import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from Assistant.Time_conversion import secs_to_hr_min_sec
from Assistant.Decorate_axes import decorate_axes_L as da
from parameters import (
    t, n, l, m,                                                              # time, initial state 
    RED, CYAN, GREEN, WHITE, YELLOW, RESET,                                  # ANSI colors
    show_E_field, print_serial_prog, confined,                               # booleans
    N, L, r_max, L_map, k_max, l_max, r0, dt, time_step,                     # parameters
    evolving_atom, eta_t, SAE_model, conf_model, total_states, state_symb    # parameters
)
from functions import (
    roots, colloc_pt, theta_k,                                                         # collocation arrays
    print_grid_info, print_smat_info, print_atom_info, print_laser_info,               # verbous printing
    f, g_lm, Y_lm, conf_selector, Absorber_func, E_field, V_int, dipole_moment, Ps     # functions
)

parser = argparse.ArgumentParser()
parser.add_argument('-v', action='store_true', help='verbose mode print')
parser.add_argument('-A', '--auto', action='store_true', help='enable automatic input GPSM and S-matrix data')
parser.add_argument('--plot', action='store_true', help='allows to plot computed results')
args = parser.parse_args()

if args.v:
    print_laser_info()
else:
    print_grid_info()
    print_smat_info()
    print_atom_info()
    print_laser_info()

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

data_dir = 'Confined_atom' if confined else 'Free_atom'

"""
Specify the data file names for the 'compatible' GPSM_states and S_matrix.

    python3 time_evolution.py [-A | --auto]

With the -A (--auto) flag, the program automatically detects the input GPSM and S-matrix files 
based on the parameters defined in parameters.py and builds upon them.

[Extremely Important Notice]:
By "compatible," the parameters used in these files must match those defined 
in parameters_and_functions.py.

Example:
If parameters_and_functions.py defines:
    L_map = 80
but the chosen files are :
    state_file = 'He_States_SAE-M1_l=0_nos=10_N=200_rmax=200_Lmap=20.dat'
    S_matrix_file = 'He_Smatrix_SAE-M1_m=0_lmax=20_kmax=50_N=200_r_max=200_L_map=20_dt=0.1.dat'
then the code will give wrong results. This is because the imported nonlinear radial mapping 
function in this script from parameters_and_functions.py:
    def f(x, Lmap=L_map):
produces a radial grid that does not match the one encoded in the data files.

In short: ensure that the parameters in the data file names are consistent 
with those in parameters_and_functions.py, otherwise the grid and data 
will be incompatible.

[For crosschecking]: running this file will show what parameter values 
parameters_and_functions.py (and this script) is currently using.
Make sure these parameters are matching with the given datafile names: `state_file' and `S_matrix_file'.
"""

if args.auto:
    if not confined:
        state_file = f'{evolving_atom}_States_{SAE_model}_l={l}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.dat'
        S_matrix_file = f'{evolving_atom}_Smatrix_{SAE_model}_m={m}_lmax={l_max}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.npy'
    else:
        state_file = f'{evolving_atom}@C60_States_{SAE_model}_{conf_model}_l={l}_nos={total_states}_N={N}_rmax={r_max}_Lmap={L_map}.dat'
        S_matrix_file = f'{evolving_atom}@C60_Smatrix_{SAE_model}_{conf_model}_m={m}_lmax={l_max}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.npy'
else:
    state_file = 'H_States_SAE-M1_l=0_nos=5_N=200_rmax=200_Lmap=80.dat'
    S_matrix_file = 'H_Smatrix_SAE-M1_m=0_lmax=20_kmax=50_N=200_rmax=200_Lmap=80_dt=0.1.npy'

state_path = this_dir / 'GPSM_states_S-matrix' / 'data_GPSM_states_S-matrix' / data_dir / state_file
state_data = np.loadtxt(state_path, skiprows=1).T

S_matrix_path = this_dir / 'GPSM_states_S-matrix' / 'data_GPSM_states_S-matrix' / data_dir / S_matrix_file
S_matrix = np.load(S_matrix_path, allow_pickle=False)          # shape: (l_max+1, N-1, N-1), dtype=complex128



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Printing system info. and show Electric field.           |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if confined:
    conf_string = conf_selector(conf_model, 0)[1]

    print(f" {WHITE}~~~~~~~~~~~~~: Conf info :~~~~~~~~~~~~~{RESET}")
    parts = conf_string.split('_')                  # Split by underscores
    confmodel = parts[0]                            # The first part (before first '_') is the confinement model
    params = [p for p in parts[1:] if '=' in p]     # Remaining parts contain key=value pairs

    print(f" {'Conf model':<10}: {confmodel}")
    for p in params:
        key, val = p.split('=')
        print(f" {key:<10}: {val}")

print(f'\n {WHITE}~~~~~~~~~~~: Time Evolution :~~~~~~~~~~{RESET}')
print(' Evolving atom            :', evolving_atom)
print(f' Evolving initial state   : (n, l, m) : ({n+l}, {l}, {m}) ~', state_symb)
print(f' θ_k[0]                   : {np.round(theta_k[0] * 180/np.pi, 4)} deg')
print(f' θ_k[-1]                  : {np.round(theta_k[-1] * 180/np.pi, 4)} deg')
print(' S_matrix file name       :', S_matrix_file)
print(' Initial state file name  :', state_file)
print(' Absorber radius (r_0)    :', r0)
print(' Total time steps         :', time_step)
print(f' Estimated time (h, m, s) : {GREEN}{secs_to_hr_min_sec(eta_t * time_step)}{RESET} ~ {GREEN}{int(eta_t * time_step)}{RESET} seconds\n')

if show_E_field and args.plot:
    fig1 = plt.figure(figsize=fig_size)
    ax1 = fig1.add_subplot(111)
    da.decorate_2d(ax1)

    E_array = E_field(t)
    ax1.plot(t, E_array)
    ax1.set_title(f'E(t) (a.u.) : max step = {time_step}', fontsize=15)
    ax1.axvline(t[time_step], linestyle='--', color='royalblue', label=f't[{time_step}]')
    ax1.set_ylim(2 * min(E_array), 2 * max(E_array))
    ax1.fill_between(t, E_array, alpha=0.2)
    ax1.legend(loc='upper right', fontsize=15, framealpha=0.5, edgecolor='w')
    plt.show(block=True)    # To show plot blocking code execution, until you manually close the plot.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   Arrays to hold evolution data                    |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r = f(colloc_pt)                                        # Nonlinear Radial collocation grid (a.u)
A_r = state_data[1:][n - 1].astype(np.complex128)       # Selecting the initial state from state_data array and storing in A_r
absorber = np.array([Absorber_func(ri) for ri in r])    # The absorbing layer.

init_glm = np.zeros((l_max + 1, N - 1), dtype=np.complex128)    # Empty zero array to hold all the time evolving partial waves.
init_glm[l - m] = A_r                           # The initial state is set as the only living partial wave. (be careful with the indexing here)
#                                               # Instead of making full 3d initial state and calculating partial waves, I can directly
#         [Initializing: partial waves]         # use the A_nl(r) as the initial partial wave.
#                                               # Details in: Section~2.3.8 `Multiple-step time evolution', see Figure 2.23(a)
glm_empty = np.empty((l_max+1, N-1), dtype=np.complex128)   # Empty gl_array to be passed in gl() function.

d_t_array = np.zeros(time_step, dtype=float)                # Dipole moment array: d(t). Doesn't include initial wavefunction's dipole moment
population_den_array = np.zeros(time_step, dtype=float)     # Array to store Population density
gs_correlation = np.zeros(time_step, dtype=float)           # Array to store the correlation function with the ground state (initial state --> A_r)
zero_psi = np.zeros((L+1, N-1), dtype=np.complex128)        # To initiate the wavefunction A(ri, θj) before each loop.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    STARTING MAIN TIME EVOLUTION                    |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p_step = 10                 # p_step = progress step. print every p_step(%) completion
checkpoints = {min(int(i * time_step / 100), time_step - 1): i for i in range(0, 101, p_step)}

start_time = time.perf_counter()            # Start measuring execution time (serial time evolution)
Y_lm_cos_theta_j = np.array([[Y_lm(l_ind+m, m, roots[j]) for j in range(L+1)] for l_ind in range(l_max+1)])    # precomputed Spherical Harmonics
for ti in range(time_step):
    psi_1 = 0 * zero_psi                    # ψ1(r, θ) = exp{-iH0(dt)/2} • ψ0(r, θ)                  # See Eq.~2.84
    psi_2 = 0 * zero_psi                    # ψ2(r, θ) = exp{-iV(r, θ, t+dt/2)(dt)/2} • ψ1(r, θ)     # See Eq.~2.86

    for j in range(L+1):                    # angular grid index j
        for l_index in range(l_max):        # summation index on 'l' of Eq.~2.85.
            psi_1[j] += np.dot(S_matrix[l_index], init_glm[l_index]) * Y_lm_cos_theta_j[l_index, j]     # calculating ψ1(r, θ)
        psi_2[j] = np.exp(-1j * V_int(r, theta_k[j], t[ti]+dt/2) * dt) * psi_1[j]                       # calculating ψ2(r, θ)

    glm_tilde = g_lm(psi_2, glm_empty)                                                   # the \tilde{g}_{\ell}(r, t) of Eq.~2.88
    for l_index in range(l_max):                                                         # Again summation index on 'l' of Eq.~2.85.
        init_glm[l_index] = np.dot(S_matrix[l_index], glm_tilde[l_index]) * absorber     # Implementing Eq.~2.89 and updating init_glm with absorbing function
    # ~~~~~~~~~~~~~~~~~~: Main time evolution algorithm ends here :~~~~~~~~~~~~~~~~~

    if print_serial_prog and ti in checkpoints:
        print(f" Evolution step {YELLOW}{ti:<5}{RESET}: {GREEN}{checkpoints[ti]:6.1f}%{RESET}")          # It will show how much the process is completed.

    d_t_array[ti] = dipole_moment(r, init_glm)                # calculating and storing the dipole moment of this instant.
    population_den_array[ti] = Ps(init_glm)                   # calculating and storing the population density.
    gs_correlation[ti] = np.abs(np.sum(A_r * init_glm[0])**2) # calculating and storing the GS correlation function.

end_time = time.perf_counter()             # ending time measurement.
wall_time = end_time - start_time          # Wall time for computing the total time evolution.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Saving dipole moment; survival probability & correlation functions |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if not confined:
    Evo_data_file = f'Evo_nopt={time_step}_{evolving_atom}({state_symb})_m={m}_{SAE_model}_L={L}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.dat'
else:
    Evo_data_file = f'Evo_nopt={time_step}_{evolving_atom}({state_symb})@C60_m={m}_{SAE_model}_{conf_model}_L={L}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.dat'

d_file_path = this_dir / 'Time_evolution_data' / f'{Evo_data_file}'

header = "t(a.u.)       E(t)(a.u.)      d(t)(a.u.)      Ps(t)       C_gs(t)"
data = np.column_stack([t[:time_step], E_field(t[:time_step]), d_t_array, population_den_array, gs_correlation])
np.savetxt(d_file_path, data, header=header, comments='', fmt='%.16e')


if wall_time > 300.0:
    print(f"\n Average wall-time per step (eta_tv)     : {GREEN}{wall_time / time_step:.5f}{RESET} seconds")
    print(f' Total wall-time for all steps (h, m, s) : {GREEN}{secs_to_hr_min_sec(wall_time)}{RESET}')
else:
    print(f"\n Average wall-time per step (eta_tv) : {GREEN}{wall_time / time_step:.5f}{RESET} seconds")
    print(f' Total wall-time for all steps       : {GREEN}{wall_time:.5f}{RESET} seconds')

print(f" evo_data_file = '{YELLOW}{Evo_data_file}{RESET}'\n")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              plotting                              |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if args.plot:
    fig2 = plt.figure(figsize=(16, 9))
    ax2 = fig2.add_subplot(221)                         # Electric Field
    ax3 = fig2.add_subplot(223)                         # dipole moment
    ax4 = fig2.add_subplot(222)                         # survival probability
    ax5 = fig2.add_subplot(224)                         # correlation function
    da.decorate_2d([ax2, ax3, ax4, ax5])

    ax3.plot(d_t_array, lw=1.5, color='deeppink', label='d(t)')
    ax2.plot([E_field(ti) for ti in t[0:time_step]], lw=1.5, color='#58C4DD', label='E(t)')
    ax4.plot(population_den_array, lw=1.5, color='orangered', label=r'P$_s$(t)')
    ax5.plot(gs_correlation, lw=1.5, color='crimson', label=r'C$_{gs}$(t)')

    ax3.set_xlabel('time steps', fontsize=15)
    ax5.set_xlabel('time steps', fontsize=15)

    ax3.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
    ax2.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
    ax4.legend(loc='upper right', fontsize=15, framealpha=0.5, edgecolor='w')
    ax5.legend(loc='upper right', fontsize=15, framealpha=0.5, edgecolor='w')
    fig2.suptitle(Evo_data_file, fontsize=13)

    fig2.subplots_adjust(
        top=0.91,
        bottom=0.075,
        left=0.04,
        right=0.985,
        hspace=0.16,
        wspace=0.125
    )

    plt.show()
