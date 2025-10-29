"""
File: vector_time_evolution.py
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
- A VECTOR implementation of time_evolution.py
- Only partial waves are evolved: glm(t+dt/2) = S(l) * glm(t)
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))   # Ensure project root (HHG-SaDAS) is in sys.path

import time
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from Assistant.Time_conversion import secs_to_hr_min_sec
from Assistant.Decorate_axes import decorate_axes_L as da
from parameters_and_functions import (
    n, l, m,                                                                                                    # initial state
    t, roots, colloc_pt, theta_k,                                                                               # arrays
    show_E_field, print_serial_prog, confined,                                                                  # booleans
    f, g_lm, Y_lm, conf_selector, Absorber_func, state_name, E_field, dipole_moment, Ps,                        # functions
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
    state_file = 'He_States_SAE-M1_l=0_nos=10_N=200_rmax=200_Lmap=20.dat'
    S_matrix_file = 'He_Smatrix_SAE-M1_m=0_lmax=20_kmax=50_N=200_r_max=200_L_map=20_dt=0.1.dat'
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

state_file = 'H_States_SAE-M1_l=0_nos=10_N=200_rmax=200_Lmap=80.dat'
state_path = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_S-matrix_data' / data_dir / state_file
state_data = np.loadtxt(state_path, skiprows=1).T

S_matrix_file = 'H_Smatrix_SAE-M1_m=0_lmax=20_kmax=50_N=200_rmax=200_Lmap=80_dt=0.1.npy'
S_matrix_path = this_dir / 'GPSM_states_S-matrix' / 'GPSM_states_S-matrix_data' / data_dir / S_matrix_file
S_matrix = np.load(S_matrix_path, allow_pickle=False)          # shape: (l_max+1, N-1, N-1), dtype=complex128



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
A_r = state_data[1:][n - 1]                             # Selecting the initial state from state_data array and storing in A_r
absorber = np.array([Absorber_func(ri) for ri in r])    # The absorbing layer.

init_glm = np.zeros((l_max + 1, N - 1), dtype=np.complex128)    # Empty zero array to hold all the time evolving partial waves.
init_glm[l - m] = A_r.astype(np.complex128)     # The initial state is set as the only living partial wave. (be careful with the indexing here)
#                                               # Instead of making full 3d initial state and calculating partial waves, I can directly
#         [Initializing: partial waves]         # use the A_nl(r) as the initial partial wave.
#                                               # Details in: Section~2.3.8 `Multiple-step time evolution', see Figure 2.23(a)
glm_empty = np.empty((l_max+1, N-1), dtype=np.complex128)   # Empty gl_array to be passed in gl() function.

d_t_array = np.zeros(time_step, dtype=float)                # Dipole moment array: d(t). Doesn't include initial wavefunction's dipole moment
population_den_array = np.zeros(time_step, dtype=float)     # Array to store Population density
zero_psi = np.zeros((L+1, N-1), dtype=np.complex128)        # To initiate the wavefunction A(ri, θj) before each loop.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    STARTING MAIN TIME EVOLUTION                    |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p_step = 10                 # p_step = progress step. print every p_step(%) completion
checkpoints = {min(int(i * time_step / 100), time_step - 1): i for i in range(0, 101, p_step)}

# Precompute spherical harmonics matrix
Y_lm_cos_theta_j = np.array([[Y_lm(l_ind + m, m, roots[j]) for j in range(L + 1)] for l_ind in range(l_max + 1)])  # (l_max+1, L+1)
Y_T = Y_lm_cos_theta_j.T   # shape: (L+1, l_max+1)

# Make sure arrays friendly for BLAS and contiguous memory
S_matrix = np.ascontiguousarray(S_matrix)          # (l_max+1, N-1, N-1)
init_glm = np.ascontiguousarray(init_glm)          # (l_max+1, N-1)
absorber = np.ascontiguousarray(absorber).astype(np.complex128)  # (N-1,)
r = np.ascontiguousarray(r)                        # (N-1,)
cos_theta = np.ascontiguousarray(roots)            # (L+1,)

# Allocate temporaries once to avoid reallocation in the loop
A = np.empty_like(init_glm)                        # (l_max+1, N-1)  holds S_matrix @ init_glm for each l
psi_1 = np.empty((L + 1, init_glm.shape[1]), dtype=np.complex128)  # (L+1, N-1)
psi_2 = np.empty_like(psi_1)                       # (L+1, N-1)
glm_tilde = glm_empty                              # reuse user's provided empty array

start_time = time.perf_counter()            # Start measuring execution time (serial time evolution)

# Main time-stepping loop (vectorized)
for ti in range(time_step):
    tmid = t[ti] + dt / 2.0

    # 1) Batched matmul: A[l, :] = S_matrix[l] @ init_glm[l, :]
    #    shapes: S_matrix (l_max+1, N-1, N-1) and init_glm[...,None] (l_max+1, N-1, 1) -> (l_max+1, N-1, 1)
    A[:] = np.squeeze(np.matmul(S_matrix, init_glm[..., None]), axis=-1)

    # 2) Construct psi_1 across all angles at once:
    #    psi_1[j, r] = sum_l Y_lm_cos_theta_j[l, j] * A[l, r]
    #    Matrix multiply: (L+1, l_max+1) @ (l_max+1, N-1) -> (L+1, N-1)
    psi_1[:] = Y_T.dot(A)

    # 3) Build potential V(theta_j, r, tmid) vectorized:
    #    V_matrix[j, r] = -E_field(tmid) * r[r] * cos(theta_j)
    #    E_field(tmid) is scalar; broadcast cos_theta[:,None] * r[None,:] -> (L+1, N-1)
    E_t = E_field(tmid)
    V_matrix = -E_t * (cos_theta[:, None] * r[None, :])   # (L+1, N-1)

    # 4) Apply exp(-i V dt) factor elementwise (vectorized over j and r)
    psi_2[:] = np.exp(-1j * V_matrix * dt) * psi_1

    # 5) Project angular -> radial partial waves (your g_lm is already vectorized)
    glm_tilde = g_lm(psi_2, glm_tilde)   # returns (l_max+1, N-1)

    # 6) Update init_glm with S_matrix * glm_tilde and apply absorber (batched matmul)
    init_glm[:] = np.squeeze(np.matmul(S_matrix, glm_tilde[..., None]), axis=-1) * absorber

    # Progress printing
    if print_serial_prog and ti in checkpoints:
        print(f"Evolution step {ti:<6}: {checkpoints[ti]:6.1f}%")

    # Dipole and survival probability
    d_t_array[ti] = dipole_moment(r, init_glm)
    population_den_array[ti] = Ps(init_glm)

end_time = time.perf_counter()             # ending time measurement.
wall_time = end_time - start_time          # Wall time for computing the total time evolution.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Saving dipole moment; survival probability & correlation functions |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if not confined:
    Evo_data_file = f'VEvo_nopt={time_step}_{evolving_atom}({state_name(n + l, l)})_m={m}_{SAE_model}_L={L}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.dat'
else:
    Evo_data_file = f'VEvo_nopt={time_step}_{evolving_atom}({state_name(n + l, l)})@C60_m={m}_{SAE_model}_{confinement_model}_L={L}_kmax={k_max}_N={N}_rmax={r_max}_Lmap={L_map}_dt={dt}.dat'

output_dir = this_dir / 'Time_evolution_data'
output_dir.mkdir(parents=True, exist_ok=True)           # Create if it doesn't exist
d_file_path = output_dir / f'{Evo_data_file}'

header = "t(a.u.)       E(t)(a.u.)      d(t)(a.u.)      Ps(t)"
data = np.column_stack([t[:time_step], E_field(t[:time_step]), d_t_array, population_den_array])
np.savetxt(d_file_path, data, header=header, comments='', fmt='%.16e')


print('\n')
print(f"evo_data_file = '{Evo_data_file}'")
print('\n')

if wall_time > 300.0:
    print(f"Average wall-time per step (eta_t)      : {wall_time / time_step:.5f} seconds")
    print(f'Total wall-time for all steps (h, m, s) : {secs_to_hr_min_sec(wall_time)}')
else:
    print(f"Average wall-time per step (eta_t) : {wall_time / time_step:.5f} seconds")
    print(f'Total wall-time for all steps      : {wall_time:.5f} seconds')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              plotting                              |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig2 = plt.figure(figsize=(16, 9))
ax2 = fig2.add_subplot(221)                         # Electric Field
ax3 = fig2.add_subplot(223)                         # dipole moment
ax4 = fig2.add_subplot(122)                         # survival probability
da.decorate_2d([ax2, ax3, ax4])

ax3.plot(d_t_array, lw=1.5, color='deeppink', label='d(t)')
ax2.plot([E_field(ti) for ti in t[0:time_step]], lw=1.5, color='#58C4DD', label='E(t)')
ax4.plot(population_den_array, lw=1.5, color='orangered', label=r'P$_s$(t)')

ax3.set_xlabel('time steps', fontsize=15)
ax4.set_xlabel('time steps', fontsize=15)

ax3.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
ax2.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='w')
ax4.legend(loc='upper right', fontsize=15, framealpha=0.5, edgecolor='w')

fig2.subplots_adjust(
    top=0.91,
    bottom=0.075,
    left=0.04,
    right=0.985,
    hspace=0.16,
    wspace=0.125
)

plt.show()
