"""
File: check_Split_operator.py
Project: HHG-SaDAS
Code Description:
    - Evolution of the initial state for one time step.
    - Normalisation and probability amplitude are checked to be okay.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

References
----------
- Appendix C -- *Derivation and Consistency of S-matrix formalism*
- Section 2.3 -- *Time Evolution of atomic wavefunction interacting with external strong-field LASER*

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
import matplotlib.pyplot as plt
from Harmonic_generation_py.parameters_and_functions import *
import Assistant.Decorate_axes.decorate_axes_L as da
start_time = time.process_time()

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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Importing files                           |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
state_file = this_dir.parent / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Free_atom' / state_file
state_data = pd.read_excel(state_file, header=None, skiprows=1).to_numpy().T

S_matrix_file = 'H_Smatrix_SAE-M1__m=0_lmax=20_kmax=50_N=200_r_max=200_L_map=80_dt=0.1.xlsx'
S_matrix_full_path = this_dir.parent / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Free_atom' / S_matrix_file
S_matrix_data = pd.read_excel(S_matrix_full_path, header=None, skiprows=1).to_numpy().T
S_matrix = np.array([[complex(*map(float, elem.split(','))) for elem in column] for column in S_matrix_data]).reshape(l_max+1, N-1, N-1)


n = 3
# Index of the state you want to check for a given l value set in parameters_and_functions.py and state_file
# This state will be considered as the initial state, for the Split Operator check here.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~: [Remember] :~~~~~~~~~~~~~~~~~~~~~~~~~~
# Principal quantum number index (n ≥ 1).
# For a given orbital angular momentum quantum number l ≥ 0,
# the actual principle quantum number = (n + l)
# Example mapping of quantum numbers (n, l, m) to standard orbitals:
# ------------------------------------------
#   l   |   n   |   m   |   Orbital / State
#   ----------------------------------------
#   0   |   1   |   0   |    1s
#   1   |   1   |   0   |    2p_{z}
#   1   |   2   |   0   |    3p_{z}
#   1   |   1   |   1   |    2p_{x}
#   1   |   2   |   1   |    3p_{x}
#   3   |   1   |   0   |    4f_{z^3}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(f'initial state       : (n, l, m) : ({n+l}, {l}, {m}) ~', state_name(n + l, l))
print('S_matrix file       :', S_matrix_file)
print('S_matrix shape      :', np.shape(S_matrix))

r = f(colloc_pt)                                      # Radial coordinate in a.u
A_r = state_data[1:][n - 1]                           # Being the eigenstate of matrix hamiltonian, we'll evolve A(r).
A_mesh, _ = np.meshgrid(A_r, theta_k)             # Initial wavefunction ~ determined by (n, l).
r_m, theta_m = np.meshgrid(r, theta_k)            # creating a meshgrid of nonlinear radial grid and angular colloc. points
U_r = A_mesh * Y_lm(l, m, np.cos(theta_m))            # U(r, θ) = A(r) • Y_lm(cosθ)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                Empty arrays to hold evolution data                 |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
U_r_recon = np.zeros((L+1, N-1), dtype=np.complex128)       # Reconstructed U(r, θ) = A(r) • P_l(cosθ)
psi_1 = np.zeros((L+1, N-1), dtype=np.complex128)           # ψ1(r, θ) = exp{-iH0(dt)/2} • ψ0(r, θ)
psi_2 = np.zeros((L+1, N-1), dtype=np.complex128)           # ψ2(r, θ) = exp{-iV(r, θ, t+dt/2)(dt)/2} • ψ1(r, θ)
psi_evolved = np.zeros((L+1, N-1), dtype=np.complex128)     # ψ(r, θ, t+dt) = exp{-iH0(dt)/2} • ψ2(r, θ)
gl_empty = np.empty((l_max+1, N-1), dtype=np.complex128)      # Empty gl_array to be passed in gl() function.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Single step evolution                        |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y_lm_cos_theta_j = np.array([[Y_lm(l_ind+m, m, roots[j]) for j in range(L+1)] for l_ind in range(l_max+1)])
gl_0_array = g_lm(U_r, gl_empty)                                                      # [STEP-0]: partial wave expansion of the initial wavefunction.
for j in range(L+1):                          # looping over all angular colloc grid
    for l_index in range(l_max+1):            # looping over all partial waves
        U_r_recon[j] += gl_0_array[l_index] * Y_lm_cos_theta_j[l_index, j]            # This will confirm whether the partial wave expansion was correctly done or not.
        psi_1[j] += G(S_matrix, gl_0_array, l_index) * Y_lm_cos_theta_j[l_index, j]   # [STEP-1]: calculating psi_1 (details: Section-2.3.7)
    psi_2[j] = np.exp(-1j * V_int(r, theta_k[j], t[0] + dt / 2) * dt) * psi_1[j]      # [STEP-2]: calculating psi_2 (details: Section-2.3.7)

glm_tilde = g_lm(psi_2, gl_empty)                                                     # again calculating partial waves after interaction term being applied.
for j in range(L+1):                          # looping over all angular colloc grid
    for l_index in range(l_max+1):            # looping over all partial waves
        psi_evolved[j] += G(S_matrix, glm_tilde, l_index) * Y_lm_cos_theta_j[l_index, j]   # [STEP-3]: final evolution (details: Section-2.3.7)



print(f'max(|A[{n-1}](t=0)|^2)  :', np.max(np.abs(U_r) ** 2))          # These two should remain same.
print(f'max(|A[{n-1}](t=dt)|^2) :', np.max(np.abs(psi_evolved)**2))    # as at first step field amplitude is very nearly equal to zero.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              plotting                              |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig1 = plt.figure(figsize=fig_size)
fig2 = plt.figure(figsize=fig_size)
fig3 = plt.figure(figsize=fig_size)
fig4 = plt.figure(figsize=fig_size)
fig6 = plt.figure(figsize=fig_size)
ax1 = fig1.add_subplot(111)             # ψ0(r, θ, t=0): Radial distribution -> A(r)
ax2 = fig2.add_subplot(211)             # g0_l
ax3 = fig2.add_subplot(223)             # dr
ax4 = fig2.add_subplot(224)             # dxj
ax5 = fig3.add_subplot(121)             # ψ0(r, θ, t=0) expanded in Legendre Polynomial
ax6 = fig3.add_subplot(122)             # ψ1(r, θ) = exp{-iH0(dt)/2} • ψ0(r, θ, t=0)
ax7 = fig4.add_subplot(121)             # ψ2(r, θ) = exp{-iV(r, θ, t+dt/2)(dt)/2} • ψ1(r, θ)
ax8 = fig4.add_subplot(122)             # ψ(r, θ, t+dt) = exp{-iH0(dt)/2} • ψ2(r, θ)
ax9 = fig6.add_subplot(111)             # g_lm(r, t=0) vs g_lm(r, dt/2)

da.decorate_polar([ax1, ax5, ax6, ax7, ax8], axis=True, grid=False)
da.decorate_2d([ax2, ax3, ax4, ax9])

for l_index in range(l_max):
    y = np.abs(gl_0_array[l_index])**2
    if np.any(y > 1e-3):  # treat as nonzero, set threshold
        ax2.plot(r, y, 'o-', label=fr'l={l_index + m}', zorder=2)         # so that the surviving partial wave can be identified.
    else:
        ax2.plot(r, y, '-', label=f'l={l_index+m}', zorder=1, alpha=0.5)  # faded + behind


"""
Convert polar coordinates (r, theta) to Cartesian for plotting.
Normally: 
    x = r*cos(theta), y = r*sin(theta) 
    -> theta=0 along the horizontal axis (x-axis), angles increase counterclockwise.
Here we want:
    - theta=0 to point along the vertical direction (z-axis) as per the axis arrangement in the thesis
    - angles to increase anticlockwise from this vertical "north" direction
Therefore, we set:
    x = -r*sin(theta)  -> horizontal axis (negative sign ensures anticlockwise rotation from vertical)
    z =  r*cos(theta)  -> vertical axis
With this convention:
    - theta = 0 points straight up along z (north)
    - theta = pi/2 points along the horizontal -x direction (anticlockwise from vertical)
[honestly this docstring is written by ChatGPT]
"""
cmap = 'nipy_spectral_r'
sc = ax1.scatter(-r_m * np.sin(theta_m), r_m * np.cos(theta_m), c=np.abs(U_r)**2, s=20, cmap=cmap)
ax5.contourf(-r_m * np.sin(theta_m), r_m * np.cos(theta_m), np.abs(U_r_recon)**2, 200, cmap=cmap)
ax6.contourf(-r_m * np.sin(theta_m), r_m * np.cos(theta_m), np.abs(psi_1)**2, 200, cmap=cmap)
ax7.contourf(-r_m * np.sin(theta_m), r_m * np.cos(theta_m), np.abs(psi_2)**2, 200, cmap=cmap)
ax8.contourf(-r_m * np.sin(theta_m), r_m * np.cos(theta_m), np.abs(psi_evolved)**2, 200, cmap=cmap)
cbar = fig1.colorbar(sc, ax=ax1)

ax3.plot(np.diff(r), 'o-', label=r'dr$_i$', color=da.mc.C_L[5])
ax4.plot(np.diff(roots), 'o-', label=r'dx$_j$', color=da.mc.C_L[2])

ax9.plot(r, np.abs(gl_0_array[l-m])**2, 'o-', markersize=10, label=rf'|g(r, {l=}, {m=}, t=0)|$^2$')
ax9.plot(r, np.abs(G(S_matrix, gl_0_array, l-m))**2, 'o--', color='m', markersize=5, label=rf'|g(r, {l=}, {m=}, t=dt/2)|$^2$')

state_name = generate_states(l)[n-1]
ax1.set_title(rf'|U(r$_i$, θ$_j$)|$^2$ = |A(r$_i$) $\cdot$ Y$_{{\ell}}^{{m}}$(cosθ$_j$)|$^2$: {evolving_atom}({state_name}, {m=}); On collocation grid [r(x$_i$), $\theta_j$]', pad=20, fontsize=15)
ax2.set_title(rf"U(r$_i$, θ$_j$) = A(r$_i$) $\cdot$ Y$_{{\ell}}^{{m}}$(cosθ$_j$): {evolving_atom}({state_name}, {m=}); Initial state's partial waves (t=0)", pad=20, fontsize=15)
ax5.set_title(rf'U(r$_i$, θ$_j$) = $\sum_{{\ell=m}}^{{\ell_{{max}}+m}}$g$_{{{{\ell}}}}$(r$_i$) N$_{{\ell m}}$ P$_{{\ell m}}$(cosθ$_j$); m={m}', pad=30, fontsize=15)
ax6.set_title(r'ψ$_1$(r$_i$, θ$_j$, t+dt/2) = exp{-iH$_0$(dt)/2} $\cdot$ U(r$_i$, θ$_j$)', pad=30, fontsize=15)
ax7.set_title(r'ψ$_2$(r$_i$, θ$_j$, t+dt/2) = exp{-iV(r$_i$, θ$_j$, t+dt/2) dt/2} $\cdot$ ψ$_1$(r$_i$, θ$_j$, t+dt/2)', pad=30, fontsize=15)
ax8.set_title(r'ψ(r$_i$, θ$_j$, t+dt) = exp{-iH$_0$(dt)/2} $\cdot$ ψ$_2$(r$_i$, θ$_j$, t+dt/2)', pad=30, fontsize=15)
ax9.set_title(r'g$_{\ell m}$(r, dt/2) = S($\ell$) $\cdot$ g$_{\ell m}$(r, t=0)', pad=30, fontsize=15)

ax1.set_xlabel('X-axis (a.u.)', fontsize=15); ax1.set_ylabel('Z-axis (a.u.)', fontsize=15)
ax5.set_xlabel('X-axis (a.u.)', fontsize=15); ax5.set_ylabel('Z-axis (a.u.)', fontsize=15)
ax6.set_xlabel('X-axis (a.u.)', fontsize=15); ax6.set_ylabel('Z-axis (a.u.)', fontsize=15)
ax7.set_xlabel('X-axis (a.u.)', fontsize=15); ax7.set_ylabel('Z-axis (a.u.)', fontsize=15)
ax8.set_xlabel('X-axis (a.u.)', fontsize=15); ax8.set_ylabel('Z-axis (a.u.)', fontsize=15)
ax2.set_xlabel(r'r(x$_i$)', fontsize=15); ax2.set_ylabel(rf'|g(r, $\ell$)|$^2$', fontsize=15)
ax9.set_xlabel('r(x)', fontsize=15); ax9.set_ylabel(rf'|g(r, $\ell$, m, t)|$^2$', fontsize=15)
ax3.set_xlabel('interval index', fontsize=15)
ax4.set_xlabel('interval index', fontsize=15)

ax2.legend(loc='upper right', ncol=3, fontsize=12, framealpha=0.5, edgecolor='w')
ax9.legend(loc='best', fontsize=12, framealpha=0.5, edgecolor='w')
ax3.legend(loc='best', fontsize=15, framealpha=0.5, edgecolor='w')
ax4.legend(loc='best', fontsize=15, framealpha=0.5, edgecolor='w')

lim = 55
ax1.axis([-lim, lim, -lim, lim])
ax5.axis([-lim, lim, -lim, lim])
ax6.axis([-lim, lim, -lim, lim])
ax7.axis([-lim, lim, -lim, lim])
ax8.axis([-lim, lim, -lim, lim])
ax2.set_xlim(-1, lim)
ax9.set_xlim(-1, lim)

fig2.subplots_adjust(
    top=0.925,
    bottom=0.075,
    left=0.08,
    right=0.975,
    hspace=0.2,
    wspace=0.225
)

fig3.subplots_adjust(
    top=0.88,
    bottom=0.11,
    left=0.075,
    right=0.975,
    hspace=0.2,
    wspace=0.28
)

fig4.subplots_adjust(
    top=0.88,
    bottom=0.11,
    left=0.075,
    right=0.975,
    hspace=0.2,
    wspace=0.28
)

end_time = time.process_time()
print(f'\nCPU Time : {end_time - start_time:.2f} seconds')
plt.show()
