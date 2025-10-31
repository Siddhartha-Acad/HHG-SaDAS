"""
File: check_GPSM.py
Project: HHG-SaDAS
Code Description:
    | ***[Demonstration purpose code]***
    | This script demonstrates the GPSM eigenstates and eigenvalues in detail.
    | All parameters are imported from `parameters_and_functions.py`.
    | Users are encouraged to run this script to verify that GPSM is functioning correctly.


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- r = f(x): nonlinear mapped radial coordinate.
- Plot the results of the Pseudospectral method.
- It doesn't generate any file.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import time
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from Atomic_units import Energy_0
from Assistant.Decorate_axes import decorate_axes_L as da
from Harmonic_generation_py.parameters import *
from Harmonic_generation_py.functions import *
start_time = time.perf_counter()

# ~~~~~~~~~~~~~~~~~~~~~~: Common Figure Settings :~~~~~~~~~~~~~~~~~~~~~
width = 6.2                         # Width in inches
height = 3                          # Height in inches
fig_scale_factor = 2                # big=2 ; medium=1.5; small=1
tickslabel_size = 18
label_fontsize = 19
fig_size = (fig_scale_factor*width, fig_scale_factor*height)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def Gauss_Lobatto_quadrature(function):
    """
    integration using Gauss–Legendre quadrature.

    :param function: Function values evaluated at the Gauss–Legendre nodes or collocation points (array-like).
    :return: Approximate integral computed using the quadrature rule (float).
    """
    return np.sum(function * int_w)

# def norm_factor(function):
#     """
#     Computes the normalization factor for a wavefunction.
#     [Note]: This function is currently unused.
#
#     :param function: Wavefunction values (array-like).
#     :return: Normalization factor (float).
#     """
#     return 1 / np.sqrt(Gauss_Lobatto_quadrature(np.abs(function)**2))
#
# def normalize(function):
#     """
#     Normalize a wavefunction.
#     [Note]: This function is currently unused.
#
#     :param function: Wavefunction values (array-like).
#     :return: Normalized wavefunction (array-like).
#     """
#     return norm_factor(function) * function



r = f(colloc_pt)                    # radial coordinate in atomic unit. (Nonlinearly discretised)
r_nm = r * a0 * 10**9               # radial coordinate in nanometer (nm)
v = V_eff(l, r)                     # Effective potential: including centrifugal term (a.u.)


# ~~~~~~~~~~~~~~~~~~~~~: H matrix, Eigenvalues (E), Eigenvectors (A) :~~~~~~~~~~~~~~~~~~~~~
H_matrix = np.zeros((N - 1, N - 1))

for i in range(N - 1):
    for j in range(i, N - 1):                   # Only computing the upper triangle
        H_matrix[i, j] = H(l, i, j, model=SAE_model)

H_matrix += np.triu(H_matrix, 1).T              # Mirror to lower triangle (excluding diagonal)

E, A = eigh(H_matrix, subset_by_index=[0, 5]); A = A.T


Ip = -E[0]                          # Ionisation potential (a.u.)
Up_au = Up(E0_au, w0)               # Ponderomotive force  (a.u.)
N_cut = N_cutoff(Ip, Up_au, w0)         # Cutoff position


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: u(r) ; φ(r) :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
u_r = np.zeros(np.shape(A)); E_n = len(E)
phi = np.zeros(np.shape(A))
norm_fact = np.sqrt(N*(N+1)/2)
for Eth in range(E_n):
    for i in range(N - 1):
        phi[Eth][i] = A[Eth][i] * P_N(colloc_pt[i])
        u_r[Eth][i] = phi[Eth][i] / np.sqrt(f_p(colloc_pt[i]))
    phi[Eth] *= norm_fact                                     # φ(r) = A(r) * P_N(r)
    u_r[Eth] *= norm_fact                                     # u(r) = r * R(r)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: S matrix :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
S_matrix = (A.T * np.exp(-1j * E * dt / 2)) @ A         # Vector computation of the S-matrix

along_S_diag = sum(phi[k]**2 for k in range(E_n))


print('~~~~~~~~~~~~~~: Spectra :~~~~~~~~~~~~~~')
print(f'Ip (a.u)              : {Ip:.3f}')
print(f'Up (a.u)              : {Up_au:.3f}')
print(f'N_cutoff              : {N_cut:.3f}')
print(f'Keldysh parameter (γ) : {Keldysh(Ip, Up_au):.3f}\n')

print('~~~~~~~~~~~~~~~~: GPSM :~~~~~~~~~~~~~~~')
print('No of Eigenvalues found  :', E_n)
print('H shape                  :', np.shape(H_matrix))
print('H[0][0]                  :', H_matrix[0][0])
print('H[-1][-1]                :', H_matrix[-1][-1])
print(f'E[0](eV) ~ l={l}           :', E[0] * Energy_0)
print(f'A[{n-1}][0]                  :', A[n-1][0])
print(f'A[{n-1}][-1]                 :', A[n-1][-1])
print(f'norm A[{n-1}] = sum(A^2)     :', sum(A[n-1]**2))
print(f'norm φ[{n-1}] = int(|φ|^2)   :', Gauss_Lobatto_quadrature(phi[n - 1] ** 2))
[print(f'E[{i}]~{state_name(i+1+l, l)}'.ljust(len(str(E_n - 1) + f'~{state_name(E_n, l)}') + 3) + f' : {E[i]:<17.15f} a.u') for i in range(E_n)]

print('u(r) dataset shape       :', np.shape(u_r))
print(f'S(l={l})-matrix shape      :', np.shape(S_matrix), '\n')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: PLOTTING :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig1, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [0.35, 0.65]}, figsize=fig_size)
fig2 = plt.figure(figsize=fig_size); fig3 = plt.figure(figsize=fig_size)
fig4 = plt.figure(figsize=fig_size); fig5 = plt.figure(figsize=fig_size); fig6 = plt.figure(figsize=fig_size)
fig7 = plt.figure(figsize=fig_size); fig8 = plt.figure(figsize=fig_size); fig9 = plt.figure(figsize=fig_size)

ax1 = axs[1]                      # Distinct states stacked
ax2 = axs[0]                      # Energy level plot
ax1_twin = ax1.twinx()            # Distinct states stacked

ax3 = fig2.add_subplot(211)       # u(r)
ax4 = fig2.add_subplot(212)       # |u(x)|^2

ax5 = fig3.add_subplot(211)       # A(x)
ax6 = fig3.add_subplot(212)       # |A(x)|^2

ax7 = fig4.add_subplot(121)       # H_matrix
ax8 = fig4.add_subplot(322)       # A(r) & P$_N$(r)
ax9 = fig4.add_subplot(324)       # A(r) * P$_N$(r)
ax10 = fig4.add_subplot(326)      # φ(r)

ax11 = fig5.add_subplot(211)      # mapping
ax12 = fig5.add_subplot(234)      # dx
ax13 = fig5.add_subplot(235)      # weights of Legendre-Gauss Quadrature of order L+1
ax14 = fig5.add_subplot(236)      # dr

ax15 = fig6.add_subplot(121)      # S_matrix real
ax16 = fig6.add_subplot(122)      # S_matrix imag

ax17 = fig7.add_subplot(121)      # S_matrix mod
ax18 = fig7.add_subplot(122)      # S_matrix mod

ax19 = fig8.add_subplot(122)      # describe S matrix behaviour
ax20 = fig8.add_subplot(121)      # potential_anim maximum deformation

ax21 = fig9.add_subplot(111)      # Electric Field

da.decorate_2d([ax1, ax2, ax3, ax5, ax6, ax4, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax19, ax20, ax21])
da.decorate_2d(ax1_twin, grid=False)
da.decorate_2d([ax7, ax15, ax17, ax18, ax16], axis_ticks=False, grid=False, visible_spine='none')
ax1.tick_params(color='red', width=5, labelbottom=True, labeltop=False, labelleft=True, labelright=False)
ax1_twin.tick_params(color='red', width=5, labelbottom=True, labeltop=False, labelleft=False, labelright=True)


states_list = generate_states(l)[0:E_n]
offset = 0; offset_arr = []; amp = 1
for Eth, E_egVal in enumerate(E):
    ax2.axhline(E_egVal, xmin=0, xmax=0.5, color='#58C4DD', label=f'l={l} levels' if Eth == 0 else '', zorder=2)

    plot_state = amp * u_r[Eth] ** 2
    if Eth != 0:
        amp = max(u_r[0] ** 2) / max(u_r[Eth] ** 2)
        plot_state = amp * u_r[Eth] ** 2
        offset -= max(plot_state) + abs(min(plot_state)) + 0.05
    offset_arr.append(offset)
    ax1.plot(r, plot_state + offset, 'o-', markersize=4, label=states_list[Eth])

    # Fill between the current curve and the previous curve or horizontal line
    if Eth == 0: ax1.fill_between(r, plot_state, alpha=0.2)
    else: ax1.fill_between(r, plot_state + offset, offset, alpha=0.2)


#       ~~~~~~~~~~~~~: State number starts from 1 (n quantum number):~~~~~~~~~~~~~
ax3.plot(r_nm, u_r[n - 1], 'o-', label=states_list[n - 1])
ax4.plot(r_nm, u_r[n - 1] ** 2, 'o-', color='#83C167', label=states_list[n - 1])

ax5.plot(colloc_pt, A[n - 1], 'o-', label=states_list[n - 1])
ax6.plot(colloc_pt, A[n - 1] ** 2, 'o-', color='#83C167', label=states_list[n - 1])

ax3.fill_between(r_nm, u_r[n - 1], alpha=0.2)
ax4.fill_between(r_nm, u_r[n - 1] ** 2, color='#83C167', alpha=0.2)

ax2.plot(r_nm, v, color='m', label='V(x)')
ax2.fill_between(r_nm, v, color='m', alpha=0.10)
ax2.axis([min(r_nm), max(r_nm), -0.55, 0.05])

# interpolation = 'spline36'
interpolation = 'none'
ax7.imshow(H_matrix, interpolation=interpolation, cmap='nipy_spectral_r')
ax7.set_xlim(-0.5, 20); ax7.set_ylim(20, -0.5)

ax8.plot(colloc_pt, A[n - 1], 'o-', markersize=4, label='A(r) ~ ' + states_list[n - 1])
ax8.plot(colloc_pt, P_N(colloc_pt), 'o-', markersize=4, label=r'P$_N$(r)')
ax9.plot(colloc_pt, A[n - 1] * P_N(colloc_pt), 'o-', markersize=4, color='#83C167', label=r'φ(r)=A(r)$\cdot$P$_N$(r) ~ ' + states_list[n - 1])
ax10.plot(colloc_pt, (norm_fact * A[n - 1] * P_N(colloc_pt)) ** 2, 'o-', markersize=4, color='m', label=r'|φ(r)|$^2$ ~ ' + states_list[n - 1])
ax9.fill_between(colloc_pt, A[n - 1] * P_N(colloc_pt), color='#83C167', alpha=0.2)
ax10.fill_between(colloc_pt, (norm_fact * A[n - 1] * P_N(colloc_pt)) ** 2, color='m', alpha=0.2)

L_map_array = np.arange(10, 120, 20)
for L_map_val in L_map_array:
    ax11.plot(colloc_pt, f(colloc_pt, Lmap=L_map_val), 'o-', label=f'L_map={L_map_val}')
    ax14.plot(np.diff(f(colloc_pt, Lmap=L_map_val)), 'o-', label='dr')

ax12.plot(np.diff(colloc_pt), 'o-', color='#83C167', label='dx')
ax13.plot(weights, 'o-', color='m', label=f'LG Quadrature weights (n={L + 1})')

ax15.imshow(np.real(S_matrix), interpolation=interpolation, cmap='jet')
ax16.imshow(np.imag(S_matrix), interpolation=interpolation, cmap='jet')
ax17.imshow(np.abs(S_matrix), interpolation=interpolation, cmap='nipy_spectral_r')
ax18.imshow(np.abs(S_matrix)**2, interpolation=interpolation, cmap='nipy_spectral_r')

ax19.plot(r_nm, along_S_diag, 'o-', label=rf'$\sum_k$  ψ$_k$(r) $\cdot$ ψ$_k$(r) ~ l={l}')
ax19.fill_between(r_nm, along_S_diag, alpha=0.2)

E_array = [E_field(ti) for ti in t]
ax21.plot(t / T, E_array, lw=1.5)
ax21.plot(t[np.argmax(E_array)] / T, max(E_array), 'ro')
ax21.set_ylim(2 * min(E_array), 2 * max(E_array))
ax21.fill_between(t / T, E_array, alpha=0.2)

max_deformed_pot =  v - E0_au * r
ax20.plot(r, max_deformed_pot, color='m')
ax20.fill_between(r, max_deformed_pot, 0, color='m', alpha=0.15)
for E_i in E:
    ax20.axhline(E_i)
ax20.axis([min(r), 30, -1, 0.1])


ax1.legend(loc='upper right', ncol=2, columnspacing=0.5, fontsize=13, framealpha=0.9, edgecolor='w')
ax2.legend(loc='lower right', fontsize=12, framealpha=0.8, edgecolor='w')
ax3.legend(loc='upper right', fontsize=15, framealpha=0.8, edgecolor='w')
ax4.legend(loc='upper right', fontsize=15, framealpha=0.8, edgecolor='w')
ax5.legend(loc='upper right', fontsize=15, framealpha=0.8, edgecolor='w')
ax6.legend(loc='upper right', fontsize=15, framealpha=0.8, edgecolor='w')
ax8.legend(loc='lower center', ncol=2, fontsize=12, framealpha=0.8, edgecolor='w')
ax9.legend(loc='lower right', fontsize=12, framealpha=0.8, edgecolor='w')
ax10.legend(loc='upper right', fontsize=12, framealpha=0.8, edgecolor='w')
ax11.legend(loc='upper left', ncol=2, fontsize=12, framealpha=0.8, edgecolor='w')
ax12.legend(loc='upper left', fontsize=12, framealpha=0.8, edgecolor='w')
ax13.legend(loc='lower center', fontsize=12, framealpha=0.8, edgecolor='w')
ax14.legend(loc='upper left', fontsize=12, framealpha=0.8, edgecolor='w')
ax19.legend(loc='upper left', fontsize=15, framealpha=0.8, edgecolor='w')

ax2.set_title(f'Potential & Energy levels; l={l}', fontsize=15)
ax1.set_title(r'Radial functions: |u(r)|$^2$ ; [u(r) = rR(r)]', fontsize=15)

ax1.set_xlabel('r (a.u) :' + r'$\longrightarrow$', fontsize=15)
ax2.set_xlabel('r (nm) :' + r'$\longrightarrow$', fontsize=15)
ax2.set_ylabel(r'E$_n$ (a.u) :' + r'$\longrightarrow$', fontsize=15)
ax3.set_ylabel('u(r) :' + r'$\longrightarrow$', fontsize=15)
ax4.set_xlabel('r (nm) :' + r'$\longrightarrow$', fontsize=15)
ax4.set_xlabel('r (nm) :' + r'$\longrightarrow$', fontsize=15)
ax4.set_ylabel(r'|u(r)|$^2$ :' + r'$\longrightarrow$', fontsize=15)
ax5.set_ylabel('A(x) :' + r'$\longrightarrow$', fontsize=15)
ax6.set_xlabel('x :' + r'$\longrightarrow$', fontsize=15)
ax6.set_ylabel(r'|A(x)|$^2$ :' + r'$\longrightarrow$', fontsize=15)
ax7.set_xlabel(f'H matrix; (l={l}); interpolation: {interpolation}', fontsize=15)
ax10.set_xlabel('x :' + r'$\longrightarrow$', fontsize=15)
ax11.set_xlabel(r'collocation points (x$_i$)' + r'$\longrightarrow$', fontsize=15)
ax11.set_ylabel(r'r$_i$=f(x$_i$)' + r'$\longrightarrow$', fontsize=15)
ax15.set_xlabel(r'$\mathfrak{Re}$(S)', fontsize=20)
ax16.set_xlabel(r'$\mathfrak{Im}$(S)', fontsize=20)
ax17.set_xlabel(f'|S({l})|', fontsize=20)
ax18.set_xlabel(rf'|S({l})|$^2$', fontsize=20)
ax21.set_xlabel(r't $\longrightarrow$', fontsize=15)
ax21.set_ylabel(r'E(t) $\longrightarrow$', fontsize=15)

fig6.suptitle(f'S_matrix; (l={l}, k_max={E_n}, dt={dt}); interpolation: {interpolation}', fontsize=20)
fig7.suptitle(f'S_matrix; (l={l}, k_max={E_n}, dt={dt}); interpolation: {interpolation}', fontsize=20)

n_values = ["n=" + str(n_ind + l) for n_ind in range(1, len(u_r) + 1)]
ax1.set_yticks(offset_arr, n_values, fontsize=12)
energy_labels = [f'{np.round(value * 27.21, 2)} eV' for value in E]

ax1_twin.set_yticks(offset_arr, energy_labels, fontsize=12)
ax1_twin.set_ylim(ax1.get_ylim())

fig1.subplots_adjust(top=0.912, bottom=0.097, left=0.061, right=0.934, wspace=0.107)
fig2.subplots_adjust(top=0.935, bottom=0.09, left=0.109, right=0.928, hspace=0.112)
fig3.subplots_adjust(top=0.935, bottom=0.09, left=0.109, right=0.928, hspace=0.112)
fig4.subplots_adjust(top=0.96, bottom=0.09, left=0.00, right=0.98, wspace=0.051, hspace=0.16)
fig5.subplots_adjust(top=0.924, bottom=0.078, left=0.069, right=0.975, wspace=0.116, hspace=0.251)
fig8.subplots_adjust(top=0.924, bottom=0.078, left=0.069, right=0.975, wspace=0.116, hspace=0.251)
fig6.subplots_adjust(left=0.02, right=0.98, wspace=0.05)
fig7.subplots_adjust(left=0.02, right=0.98, wspace=0.05)

end_time = time.perf_counter()
print(f'Execution Wall-Time : {end_time - start_time:.2f} seconds')

plt.show()
