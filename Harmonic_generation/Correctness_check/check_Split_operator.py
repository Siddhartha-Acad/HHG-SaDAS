"""
~ Evolve_dt.py

• r = f(x): nonlinear mapped radial coordinate.
• Evolution of the initial state for one time step.
• Normalisation, probability amplitude are checked to be okay.
"""

import time
import pandas as pd
from scipy.special import factorial, lpmv
import matplotlib.pyplot as plt
from Harmonic_generation.parameters_and_functions import *
import Assistant.Decorate_axes.decorate_axes_L as da
start_time = time.time()

def H_psi(R, theta):
    """
    Hydrogen atom wavefunctions; will be used as the initial state (t=0).
    ψ(r, θ) = R(r).P(cos(θ))
    R: Radial function.
    θ: theta, uniform theta grid distribution.
    """
    return R * Y_lm(l, m, np.cos(theta))

def C_fact(l, m):
    """
    Orthogonality constant factor of P_lm
    """
    return 2 * factorial(l+m) / ((2*l+1)*factorial(l-m))

def N_fact(l, m):
    """
    Normalization constant of Y_lm
    """
    return (-1)**m * np.sqrt((2*l+1) * factorial(l-m) / (4*np.pi * factorial(l+m)))

def a_legendre(l, m, x):
    """
    Associated Legendre polynomial without Condon-Shortley phase
    """
    return lpmv(m, l, x) * (-1)**m

def Y_lm(l, m, x):
    """
    Y_lm(x) = N_lm * P_lm(x)
    """
    return N_fact(l, m) * a_legendre(l, m, x)

def gl(Psi_t):
    """
    Calculate partial waves for a given wavefunction.
    Psi_t : Complex array representing the wavefunction at each combination of (r, θ, t).
    Notes : θ must include θ_k, as they are required for the calculation.
    ~~ Psi_t[theta_k_locs[k]][i] = ψ(r_i, θ_k, t)
    """
    gl_arr = []
    psi_vals = np.array([[Psi_t[k][i] for k in range(len_k)] for i in range(len_r)])
    for l_ind in range(L):
        g_l = (weights * a_legendre_vals[l_ind, :]) @ psi_vals.T
        gl_arr.append(g_l / (N_fact(l_ind + m, m) * C_fact(l_ind + m, m)))
    return gl_arr

def G(Sl_matrix, gl_arr, l_ind):
    """
    G(l) = S(l) * g(l)
    """
    return np.dot(Sl_matrix[l_ind], gl_arr[l_ind])


m = 1; l = 1; n = 2
print(f'Evolving initial state   : (n, l, m) : ({n+l}, {l}, {m}) ~', state_name(n + l, l))

# ~~~~~~~~~~~~~~~~~~~~~~~: Importing files :~~~~~~~~~~~~~~~~~~~~~~~
this_dir = Path(__file__).resolve().parent


psi_file = 'He_States__l=0_nos=10_N=200_rmax=200_Lmap=20.xlsx'
psi_file = this_dir.parent / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Free_atom' / psi_file
psi_data = pd.read_excel(psi_file, header=None, skiprows=1).to_numpy().T

S_matrix_file = 'He_Smatrix__m=0_lmax=20_kmax=50_N=200_r_max=200_L_map=20_dt=0.1.xlsx'
S_matrix_file = this_dir.parent / 'GPSM_states_S-matrix' / 'GPSM_states_and_Smatrix_data' / 'Free_atom' / S_matrix_file
S_matrix_data = pd.read_excel(S_matrix_file, header=None, skiprows=1).to_numpy().T
S_matrix = np.array([[complex(*map(float, elem.split(','))) for elem in column] for column in S_matrix_data]).reshape(L+1, N-1, N-1)


print('S_matrix file name       :', S_matrix_file)
print('S_matrix shape           :', np.shape(S_matrix))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


r = f(colloc_pt)                                      # Radial coordinate in a.u
A_r = psi_data[1:][n-1]                               # Being the eigenstate of matrix hamiltonian, we'll evolve A(r).
R_m, _ = np.meshgrid(A_r, theta_k)                  # Initial wavefunction ~ determined by (n, l).
r_m, theta_m = np.meshgrid(r, theta_k)
psi_0 = H_psi(R_m, theta_m)                           # ψ0(r, θ) = A(r) • P_l(cosθ)
a_legendre_vals = np.array([[a_legendre(l_index+m, m, root) for root in roots] for l_index in range(L)])     # will be used in gl(r) function

len_r = len(r); len_k = len(weights)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: G(l) = S(l) * g(l) :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psi_0_recon = np.zeros((len_k, len_r), dtype=np.complex128)     # Reconstructed ψ0(r, θ) = A(r) • P_l(cosθ)
psi_1 = np.zeros((len_k, len_r), dtype=np.complex128)           # ψ1(r, θ) = exp{-iH0(dt)/2} • ψ0(r, θ)
psi_2 = np.zeros((len_k, len_r), dtype=np.complex128)           # ψ2(r, θ) = exp{-iV(r, θ, t+dt/2)(dt)/2} • ψ1(r, θ)
psi_evolved = np.zeros((len_k, len_r), dtype=np.complex128)     # ψ(r, θ, t+dt) = exp{-iH0(dt)/2} • ψ2(r, θ)
cos_theta = np.cos(theta_k); len_theta = len_k
Y_lm_cos_theta_j = np.array([[Y_lm(l_ind+m, m, cos_theta[j]) for j in range(len_k)] for l_ind in range(L)])


gl_0_array = gl(psi_0)
for j in range(len_k):
    for l_index in range(L):
        psi_0_recon[j] += gl_0_array[l_index] * Y_lm_cos_theta_j[l_index, j]
        psi_1[j] += G(S_matrix, gl_0_array, l_index) * Y_lm_cos_theta_j[l_index, j]
    psi_2[j] = np.exp(-1j * V_int(r, theta_k[j], t[0] + dt / 2) * dt) * psi_1[j]

gl_2_array = gl(psi_2)
for j in range(len_k):
    for l_index in range(L):
        psi_evolved[j] += G(S_matrix, gl_2_array, l_index) * Y_lm_cos_theta_j[l_index, j]

print(f'max(|A[{n-1}](t=0)|^2)       :', np.max(np.abs(psi_0)**2))
print(f'max(|A[{n-1}](t=dt)|^2)      :', np.max(np.abs(psi_evolved)**2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: plotting :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
fig4 = plt.figure()
fig5 = plt.figure()
ax1 = fig1.add_subplot(111)             # ψ0(r, θ, t=0): Radial distribution -> A(r)
ax2 = fig2.add_subplot(211)             # g0_l
ax3 = fig2.add_subplot(223)             # dr
ax4 = fig2.add_subplot(224)             # dθ
ax5 = fig3.add_subplot(121)             # ψ0(r, θ, t=0) expanded in Legendre Polynomial
ax6 = fig3.add_subplot(122)            # ψ1(r, θ) = exp{-iH0(dt)/2} • ψ0(r, θ, t=0)
ax7 = fig4.add_subplot(121)            # ψ2(r, θ) = exp{-iV(r, θ, t+dt/2)(dt)/2} • ψ1(r, θ)
ax8 = fig4.add_subplot(122)            # ψ(r, θ, t+dt) = exp{-iH0(dt)/2} • ψ2(r, θ)
ax9 = fig5.add_subplot(111)           # |A(r, t=0)|^2, |A(r, t=dt)|^2

da.decorate_polar([ax1, ax5, ax6, ax7, ax8])
da.decorate_2d([ax2, ax3, ax4, ax9])

for l_index in range(L):
    ax2.plot(r, gl_0_array[l_index], 'o-', label=rf'g$_{{0{l_index}}}$(r)')

ax1.contourf(r_m * np.sin(theta_m), r_m * np.cos(theta_m), psi_0**2, 150, cmap='jet')
ax5.contourf(r_m * np.sin(theta_m), r_m * np.cos(theta_m), psi_0_recon**2, 200, cmap='jet')
ax6.contourf(r_m * np.sin(theta_m), r_m * np.cos(theta_m), psi_1**2, 200, cmap='jet')
ax7.contourf(r_m * np.sin(theta_m), r_m * np.cos(theta_m), psi_2**2, 200, cmap='jet')
ax8.contourf(r_m * np.sin(theta_m), r_m * np.cos(theta_m), psi_evolved**2, 200, cmap='jet')

ax3.plot(np.diff(r), 'o-', label='dr')
ax4.plot(np.diff(theta_k), 'o-', color='m', label='dθ')

ax9.plot(r, np.abs(psi_0[0]) ** 2, 'o-', markersize=10, label=r'|A(θ[0], t=0)|$^2$')
ax9.plot(r, np.abs(psi_evolved[0]) ** 2, 'o--', color='m', label=r'|A(θ[0], t=dt)|$^2$')

ax2.legend(loc='upper right', ncol=3, fontsize=12, framealpha=0.5, edgecolor='k')
ax9.legend(loc='best', fontsize=12, framealpha=0.5, edgecolor='k')
ax3.legend(loc='upper left', fontsize=15, framealpha=0.5, edgecolor='k')
ax4.legend(loc='lower right', fontsize=15, framealpha=0.5, edgecolor='k')

state_name = generate_states(l)[n-1]
ax1.set_title(f'ψ$_0$(r, θ) ~ {state_name} ({m=})', pad=20, fontsize=15)
ax2.set_title(f'ψ$_0$(r, θ) ~ {state_name} ({m=}); Grid discretization', pad=20, fontsize=15)
ax5.set_title(fr'ψ$_0$(r$_i$, θ$_j$) = $\sum_{{\ell=m}}^{{L+m}}$g$_{{{{\ell}}}}$(r$_i$) N$_{{\ell m}}$ P$_{{\ell m}}$(cosθ$_j$); m={m}', pad=30, fontsize=15)
ax6.set_title(r'ψ$_1$(r, θ) = exp{-iH$_0$(dt)/2} • ψ$_0$(r, θ)', pad=30, fontsize=15)
ax7.set_title(r'ψ$_2$(r, θ) = exp{-iV(r, θ, t+dt/2)(dt)/2} • ψ$_1$(r, θ)', pad=30, fontsize=15)
ax8.set_title(r'ψ(r, θ, t+dt) = exp{-iH$_0$(dt)/2} • ψ$_2$(r, θ)', pad=30, fontsize=15)

lim = 10
ax1.axis([-lim, lim, -lim, lim])
ax5.axis([-lim, lim, -lim, lim])
ax6.axis([-lim, lim, -lim, lim])
ax7.axis([-lim, lim, -lim, lim])
ax8.axis([-lim, lim, -lim, lim])

ax2.set_xlim(-1, lim)
ax9.set_xlim(-1, lim)

fig1.subplots_adjust(top=0.876, bottom=0.048)
fig3.subplots_adjust(top=0.876, bottom=0.048)
fig4.subplots_adjust(top=0.876, bottom=0.048)
end_time = time.time()
print(f'Execution Time           : {end_time - start_time:.2f} seconds')

plt.show()
