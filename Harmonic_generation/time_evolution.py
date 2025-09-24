"""
~ time_evolution.py

• Initial state: m >= 0 states.
• Evolves the initial state for n time steps, controlled by `time_step`.
• New Algo: Evolve only the partial waves: gl(t+dt) = S(l) * gl(t)
            Partial wave expansion in terms of Y_lm.
            Hence, formula for gl also modified.
• Shows Electric field and dipole moment.
• Saves Electric field and dipole moment data in Excel file.
• Does not show final wavefunction after n steps.
• Calculates the evolution of the wavefunction only for theta = theta_k
• This increases efficiency and CPU time.

"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Assistant.Time_conversion import secs_to_hr_min_sec
from Assistant.Decorate_axes import decorate_axes_D as da
from Harmonic_generation.parameters_and_functions import t, roots, colloc_pt                                              # imported arrays
from Harmonic_generation.parameters_and_functions import n, l, m, cpp, E0_au, w0                                          # imported params
from Harmonic_generation.parameters_and_functions import N, l_max, r0, dt, evolving_atom                                  # imported params
from Harmonic_generation.parameters_and_functions import f, factorial, g_lm, Y_lm, Absorber_func, state_name  # imported funcns


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
da.decorate_2d(ax1)

def E_field(t):           # Single color
    return E0_au * np.sin(w0*t) * (np.sin(w0*t / (2*cpp))) ** 2

# def E_field(t):             # Double color
#     return (E0_au * np.sin(w0 * t) + E0_au_2 * np.sin(w0_2 * t)) * (np.sin(w0 * t / (2 * cpp))) ** 2

def V_int(r, theta, t):
    return -E_field(t) * r * np.cos(theta)

def C_fact(l, m):
    """
    Orthogonality constant factor of P_lm
    """
    return 2 * factorial(l+m, exact=False) / ((2*l+1)*factorial(l-m, exact=False))

def N_fact(l, m):
    """
    Normalization constant of Y_lm
    """
    return (-1)**m * np.sqrt((2*l+1) * factorial(l-m, exact=False) / (4*np.pi * factorial(l+m, exact=False)))

def alpha(l, m):
    """
    General factor for dipole moment
    """
    return 2*np.pi * N_fact(l+1, m) * N_fact(l, m) * C_fact(l, m) * (l + m + 1) / (2*l + 3)


def dipole_moment_gl(gl_arr):
    """
    gl_arr: partial waves of A(r, θ, t)
    """
    integrals = np.array([np.sum(r * np.real(np.conj(gl_arr[l_ind]) * gl_arr[l_ind + 1])) for l_ind in l_ind_arr])
    return 2 * np.sum(l_factors_d * integrals)

def Ps(gl_arr):
    """
    Survival probability
    """
    return np.sum(np.sum(np.abs(gl_arr)**2, axis=1))


eta_t = 0.03; time_step = 100

# ~~~~~~~~~~~~~~~~~~~~~~~: Importing files :~~~~~~~~~~~~~~~~~~~~~~~
S_matrix_file = 'He_Smatrix_SAE-M1__m=1_lmax=20_kmax=50_N=200_r_max=200_L_map=80_dt=0.1.xlsx'
file_S_matrix = rf'E:\Python_programs\HHG-SaDAS\Harmonic_generation\GPSM_states_S-matrix\GPSM_states_and_Smatrix_data\Free_atom\{S_matrix_file}'
S_matrix_data = pd.read_excel(file_S_matrix, header=None, skiprows=1).to_numpy().T
S_matrix = np.array([[complex(*map(float, elem.split(','))) for elem in column] for column in S_matrix_data]).reshape(l_max+1, N-1, N-1)

psi_file = 'He_States_SAE-M1__l=1_nos=10_N=200_rmax=200_Lmap=80.xlsx'
file_psi = rf'E:\Python_programs\HHG-SaDAS\Harmonic_generation\GPSM_states_S-matrix\GPSM_states_and_Smatrix_data\Free_atom\{psi_file}'
psi_data = pd.read_excel(file_psi, header=None, skiprows=1).to_numpy().T

# ~~~~~~~: Some pre-computed arrays to make calculations faster :~~~~~~~
r = f(colloc_pt)                                      # Radial coordinate in a.u
A_r = psi_data[1:][n-1]                               # Being the eigenstate of matrix hamiltonian, we'll evolve A(r).
absorber = np.array([Absorber_func(ri) for ri in r])

theta_k = np.arccos(roots)
l_ind_arr = np.arange(l_max)
l_factors_d = alpha(l_ind_arr+m, m)
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
E_array = [E_field(ti) for ti in t]
ax1.plot(t, E_array)
ax1.set_title(f'E(t) ~ max step = {time_step}', fontsize=15)
ax1.axvline(t[time_step], color='yellow', label=f't[{time_step}]')
ax1.set_ylim(2*min(E_array), 2*max(E_array))
ax1.fill_between(t, E_array, alpha=0.2)
ax1.legend(loc='upper right', fontsize=15, framealpha=0.5, edgecolor='k')
plt.show()


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


    print(f'Evolution step {ti}    : {((ti+1)/time_step)*100:.4f}%')

    dipole_mom = dipole_moment_gl(init_gl)
    population_den_array = np.append(population_den_array, Ps(init_gl))
    d_t_array = np.append(d_t_array, dipole_mom)

end_time = time.time()


# ~~~~~~~~~~~~~~~~~~~: Saving dipole moment :~~~~~~~~~~~~~~~~~~~
# In the dipole moment files where cpp is not mentioned, assumed to be cpp=60
# dipole_file_name = f'd_len_Ps_{time_step}_{evolving_atom}_{state_name(n+l, l)}_m={m}__L={L}_k_max={k_max}_N={N}_r_max={r_max}_L_map={L_map}_dt={dt}_wAb_r0={r0}_all_ok.xlsx'
# dipole_moment_data['d(t)'] = d_t_array
# dipole_moment_data['Ps(t)'] = population_den_array
# df_dipole_moment_data = pd.DataFrame(dipole_moment_data)
# d_file_path = rf'E:\Python_programs\HHG\GPSM\GPSM_Y_lm\free_SAE\Important_files\{dipole_file_name}'
# df_dipole_moment_data.to_excel(d_file_path, index=False)
# print(f"'{dipole_file_name}'")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: plotting :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig3 = plt.figure(figsize=(18, 9))
ax8 = fig3.add_subplot(311)                         # Electric Field
ax7 = fig3.add_subplot(312)                         # dipole moment, from t=0
ax9 = fig3.add_subplot(313)                         # dipole moment, from t=0
da.decorate_2d([ax7, ax8, ax9])

ax7.plot(d_t_array, lw=2, color='deeppink', label='d(t)')
ax8.plot([E_field(ti) for ti in t[0:time_step]], lw=2, color='#58C4DD', label='E(t)')
ax9.plot(population_den_array, lw=2, color='orangered', label=r'P$_s$(t)')

ax7.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='k')
ax8.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='k')
ax9.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='k')
fig3.subplots_adjust(top=0.92, bottom=0.06, right=0.97, left=0.048, hspace=0.14)
# fig3.suptitle(dipole_file_name, fontsize=13)

print('Average time  for each step     :', (end_time-start_time)/time_step, ' sec')
print('Total Execution Time (h, m, s)  :', secs_to_hr_min_sec(end_time - start_time))

plt.show()
