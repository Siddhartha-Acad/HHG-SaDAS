"""
~ H_ij_positive.py: H atom interacting with external Laser field.

• Understood the necessity of positive energy solutions of H_matrix.
• r = f(x): nonlinear mapped radial coordinate.
• Plot the results of the Pseudospectral method.
• All the variables are exported from here which are used by other programme files.
• It doesn't generate any file.

"""
import time
import warnings
from scipy.linalg import eigh
from Atomic_units import a0, Energy_0
from Harmonic_generation.function_bank import *
from Assistant.Time_conversion import secs_to_hr_min_sec
start_time = time.time()

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

print('*** RuntimeWarning     : Blocked from H_ij_positive.py ***')
print('*** DeprecationWarning : Blocked from H_ij_positive.py ***\n')

# ~~~~~~~~~~~~~~: Common Figure Settings :~~~~~~~~~~~~~~
width = 6.2                         # Width in inches
height = 3                          # Height in inches
fig_scale_factor = 2              # big=2 ; medium=1.5; small=1
tickslabel_size = 18
label_fontsize = 19
fig_size = (fig_scale_factor*width, fig_scale_factor*height)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def mapped_integration(function):
    return np.sum(function * int_w)

# def norm_factor(function):
#     return 1 / np.sqrt(mapped_integration(np.abs(function)**2))
#
# def normalize(function):
#     return norm_factor(function) * function



int_w = 2 / (N * (N + 1) * (legendre(N)(colloc_pt))**2)           # integration weights; used in mapped_integration.

theta_org = np.linspace(0, 2*np.pi, 200)                # original theta.
roots, weights = np.polynomial.legendre.leggauss(L+1)             # Gauss-Legendre Quadrature.
theta_k = np.arccos(roots)                                        # Gauss-Legendre angular collocation points
insert_positions = np.searchsorted(theta_org, theta_k)            # To find the proper position to insert theta_k
theta = np.insert(theta_org, insert_positions, theta_k)
theta_k_locs = np.where(np.isin(theta, theta_k))[0]



r = f(colloc_pt)                                # converting length unit from a.u to meter
r_nm = r * a0 * 10**9                           # radial coordinate in nanometer (nm)
v = V_eff(l, r)
n0_array = np.arange(1, L + 1 + 1, 1)           # Set the max energy search one level above the required level.
l0_levels = -1 / (2 * n0_array**2)              # Calculate energy levels up to just below the target state's energy.


# ------------------------------------: H matrix, E, A :------------------------------------
H_matrix = np.zeros((N - 1, N - 1))
for i in range(N - 1):
    for j in range(i, N - 1):                           # Only computing the upper triangle
        H_matrix[i, j] = H_matrix[j, i] = H(l, i, j)
E, A = eigh(H_matrix, subset_by_index=[0, 5]); A = A.T


Ip = 0.5               # Hard setting Ip=0.5; for Hydrogen atom ~ E0(1s). Generally, Ip = -E[0]
Up_au = Up(E0_au, w0)
N_cut = N_cutoff(Ip, Up_au)


if __name__ == '__main__':
    # ----------------------------:  u(r) = psi(f(x)) ; φ(r) :----------------------------
    psi = np.zeros(np.shape(A)); E_n = len(E)
    phi = np.zeros(np.shape(A))
    norm_fact = np.sqrt(N*(N+1)/2)
    for Eth in range(E_n):
        for i in range(N - 1):
            phi[Eth][i] = A[Eth][i] * P_N(colloc_pt[i])
            psi[Eth][i] = phi[Eth][i] / np.sqrt(f_p(colloc_pt[i]))
        phi[Eth] *= norm_fact             # φ(r) = A(r) * P_N(r)
        psi[Eth] *= norm_fact             # psi(r) = u(r) = r * R(r)

    # ------------------------------------: S matrix :------------------------------------
    S_matrix_real = np.zeros((N-1, N-1))
    S_matrix_imag = np.zeros((N-1, N-1))
    for i in range(N - 1):
        for j in range(i, N - 1):                               # Only compute the upper triangle
            matrix_ele = S(E, A, i, j)
            S_matrix_real[i][j] = np.real(matrix_ele)
            S_matrix_imag[i][j] = np.imag(matrix_ele)
            if i != j:                                          # Mirror the values to the lower triangle
                S_matrix_real[j][i] = np.real(matrix_ele)
                S_matrix_imag[j][i] = np.imag(matrix_ele)
    along_S_diag = sum(phi[k]**2 for k in range(E_n))

    print('~~~~~~~~~~~~: Spectra :~~~~~~~~~~~~')
    print(f'Ip (a.u)              : {Ip:.3f}')
    print(f'Up (a.u)              : {Up_au:.3f}')
    print(f'N_cutoff              : {N_cut:.3f}')
    print(f'Keldysh parameter (γ) : {Keldysh(Ip, Up_au):.3f}\n')

    print('~~~~~~~~~~~~: GPSM :~~~~~~~~~~~~')
    print('No of Eigenvalues found  :', E_n)
    print('H shape                  :', np.shape(H_matrix))
    print('H[0][0]                  :', H_matrix[0][0])
    print('H[-1][-1]                :', H_matrix[-1][-1])
    print(f'E[0](eV) ~ l={l}           :', E[0] * Energy_0)
    print(f'A[{n-1}][0]                  :', A[n-1][0])
    print(f'A[{n-1}][-1]                 :', A[n-1][-1])
    print(f'norm A[{n-1}] = sum(A^2)     :', sum(A[n-1]**2))
    print(f'norm φ[{n-1}] = int(|φ|^2)   :', mapped_integration(phi[n-1]**2))

    [print(f'E[{i}]~{state_name(i+1+l, l)}'.ljust(len(str(E_n - 1) + f'~{state_name(E_n, l)}') + 3) +
           f' : {E[i]:<17.15f} a.u, Rel_Error: {rel_E_error(l + i + 1, E[i]):.5e}') for i in range(E_n)]

    print('ψ shape                  :', np.shape(psi))
    print('S shape                  :', np.shape(S_matrix_real), '\n')


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: PLOTTING :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    from Assistant.Decorate_axes import decorate_axes_L as da
    import matplotlib.pyplot as plt
    fig2, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [0.35, 0.65]}, figsize=fig_size)
    fig1 = plt.figure(figsize=fig_size); fig3 = plt.figure(figsize=fig_size); fig3a = plt.figure(figsize=fig_size)
    fig4 = plt.figure(figsize=fig_size); fig5 = plt.figure(figsize=fig_size); fig6 = plt.figure(figsize=fig_size)
    fig6a = plt.figure(figsize=fig_size); fig7 = plt.figure(figsize=fig_size); fig8 = plt.figure(figsize=fig_size)
    ax1 = fig1.add_subplot(221)         # Eigenstates with potential_anim function
    ax2 = fig1.add_subplot(222)         # modSq Eigenstates with potential_anim function
    ax3 = fig1.add_subplot(234)         # Energy levels with potential_anim
    ax4 = fig1.add_subplot(235)         # Energy levels o- plots compared with E0/n^2
    ax5 = fig1.add_subplot(236)         # Energy levels bar plot
    ax8 = fig3.add_subplot(211)         # u(r)
    ax9 = fig3.add_subplot(212)         # |u(x)|^2
    ax8a = fig3a.add_subplot(211)       # A(x)
    ax8b = fig3a.add_subplot(212)       # |A(x)|^2
    ax6 = axs[1]                        # Distinct states stacked
    ax6_twin = ax6.twinx()              # Distinct states stacked
    ax7 = axs[0]                        # Energy level plot
    ax10 = fig4.add_subplot(121)        # H_matrix
    ax11 = fig4.add_subplot(322)        # A(r) & P$_N$(r)
    ax12 = fig4.add_subplot(324)        # A(r) * P$_N$(r)
    ax13 = fig4.add_subplot(326)        # φ(r)
    ax14 = fig5.add_subplot(211)        # mapping
    ax15 = fig5.add_subplot(234)        # dx
    ax15_1 = fig5.add_subplot(235)      # weights of Legendre-Gauss Quadrature of order L+1
    ax16 = fig5.add_subplot(236)        # dr
    ax17 = fig6.add_subplot(121)        # S_matrix real
    ax18 = fig6.add_subplot(122)        # S_matrix imag
    ax17_1 = fig6a.add_subplot(121)     # S_matrix mod
    ax17_2 = fig6a.add_subplot(122)     # S_matrix mod
    ax19 = fig7.add_subplot(122)        # describe S matrix behaviour
    ax21 = fig7.add_subplot(121)        # potential_anim maximum deformation
    ax20 = fig8.add_subplot(111)        # Electric Field

    da.decorate_2d([ax1, ax2, ax3, ax4, ax5, ax7, ax8, ax8a, ax8b, ax9, ax11, ax12, ax13, ax14, ax15, ax15_1, ax16, ax19, ax20, ax21])
    da.decorate_2d(ax6, tick_param=False)
    da.decorate_2d(ax6_twin, grid=False, tick_param=False)
    da.decorate_2d([ax10, ax17, ax17_1, ax17_2, ax18], axis_ticks=False, grid=False, visible_spine='none')
    ax6.tick_params(color='red', width=5, labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax6_twin.tick_params(color='red', width=5, labelbottom=True, labeltop=False, labelleft=False, labelright=True)


    states_list = generate_states(l)[0:E_n]
    offset = 0; offset_arr = []; amp = 1
    for i in range(len(l0_levels)):
        ax3.axhline(l0_levels[i], color='w', linestyle='dashed', lw=1, alpha=0.5, label='l=0 levels' if i == 0 else '')
        ax7.axhline(l0_levels[i], color='w', linestyle='dashed', lw=1, alpha=0.5, label='l=0 levels' if i == 0 else '')

    for Eth, E_egVal in enumerate(E):
        ax3.axhline(E_egVal, xmin=0, xmax=0.5, color='#58C4DD', label='E levels' if Eth == 0 else '')
        ax7.axhline(E_egVal, xmin=0, xmax=0.5, color='#58C4DD', label=f'l={l} levels' if Eth == 0 else '', zorder=2)

        ax1.plot(r, psi[Eth], label=states_list[Eth])
        ax2.plot(r_nm, abs(psi[Eth]) ** 2, label=states_list[Eth])
        # print(f'E({Eth}) = {np.round(E_egVal * 27.21, 7)} eV')
        plot_state = amp * psi[Eth] ** 2
        if Eth != 0:
            amp = max(psi[0] ** 2) / max(psi[Eth] ** 2)
            plot_state = amp * psi[Eth] ** 2
            offset -= max(plot_state) + abs(min(plot_state)) + 0.05
        offset_arr.append(offset)
        ax6.plot(r, plot_state + offset, 'o-', markersize=4, label=states_list[Eth])

        # Fill between the current curve and the previous curve or horizontal line
        if Eth == 0: ax6.fill_between(r, plot_state, alpha=0.2)
        else: ax6.fill_between(r, plot_state + offset, offset, alpha=0.2)

    #       ~~~~~~~~~~~~~: State number starts from 1 (n quantum number):~~~~~~~~~~~~~
    ax8.plot(r_nm, psi[n - 1], 'o-', label=states_list[n - 1])
    ax9.plot(r_nm, psi[n - 1] ** 2, 'o-', color='#83C167', label=states_list[n - 1])

    ax8a.plot(colloc_pt, A[n - 1], 'o-', label=states_list[n - 1])
    ax8b.plot(colloc_pt, A[n - 1] ** 2, 'o-', color='#83C167', label=states_list[n - 1])

    ax8.fill_between(r_nm, psi[n - 1], alpha=0.2)
    ax9.fill_between(r_nm, psi[n - 1] ** 2, color='#83C167', alpha=0.2)



    ax3.plot(r_nm, v, color='m', label='V(x)')
    ax1.plot(r, v, color='m')
    ax2.plot(r_nm, v, color='m')
    ax7.plot(r_nm, v, color='m', label='V(x)')
    ax1.fill_between(r, v, color='m', alpha=0.15)
    ax7.fill_between(r_nm, v, color='m', alpha=0.10)
    ax3.fill_between(r_nm, v, color='m', alpha=0.15)
    ax2.fill_between(r_nm, v, color='m', alpha=0.15)
    ax7.axis([min(r_nm), max(r_nm), -0.55, 0.05])


    n_array = np.arange(1, 80, 1)
    ref_n_array = n_array[0 : len(l0_levels)-2]
    plot_n_array = n_array[l: len(E)+l]
    ax4.plot(plot_n_array, E, 'o--', color='yellow', label=f'num (l={l})', zorder=2)
    ax4.plot(ref_n_array, l0_levels[0] / ref_n_array ** 2, 'o-', color='m', label='E' + r'$_0$' + '/n' + '$^2$', zorder=1)
    ax4.plot(ref_n_array, l0_levels[0:-2], 'o-', label='l=0 levels', zorder=0)
    ax4.set_xlim(0, max(ref_n_array)+1)

    ax1.set_ylim(-0.42, 0.78)
    ax2.set_ylim(-0.05, 0.55)
    # ax8.set_xlim(-1, 10)
    # ax9.set_xlim(-1, 10)

    ax5.bar(plot_n_array, E, color=da.dec_color[l: len(E)+l], label=f'num (l={l})')
    ax5.bar(ref_n_array, l0_levels[0:-2], color=da.dec_color[:len(l0_levels)-2], alpha=0.4, label='l=0 levels')

    # interpolation = 'spline36'
    interpolation = 'none'
    ax10.imshow(H_matrix, interpolation=interpolation, cmap='nipy_spectral_r')
    ax10.set_xlim(-0.5, 20)
    ax10.set_ylim(20, -0.5)
    # ax10.axis('off')

    ax11.plot(colloc_pt, A[n - 1], 'o-', markersize=4, label='A(r) ~ ' + states_list[n - 1])
    ax11.plot(colloc_pt, P_N(colloc_pt), 'o-', markersize=4, label=r'P$_N$(r)')
    ax12.plot(colloc_pt, A[n - 1] * P_N(colloc_pt), 'o-', markersize=4, color='#83C167', label=r'φ(r)=A(r)$\cdot$P$_N$(r) ~ ' + states_list[n - 1])
    ax13.plot(colloc_pt, N * (N + 1) / 2 * (A[n - 1] * P_N(colloc_pt)) ** 2, 'o-', markersize=4, color='m', label=r'|φ(r)|$^2$ ~ ' + states_list[n - 1])
    ax12.fill_between(colloc_pt, A[n - 1] * P_N(colloc_pt), color='#83C167', alpha=0.2)
    ax13.fill_between(colloc_pt, N * (N + 1) / 2 * (A[n - 1] * P_N(colloc_pt)) ** 2, color='m', alpha=0.2)

    L_map_array = np.arange(10, 120, 20)
    for L_map_val in L_map_array:
        ax14.plot(colloc_pt, f(colloc_pt, Lmap=L_map_val), 'o-', label=f'L_map={L_map_val}')
        ax16.plot(np.diff(f(colloc_pt, Lmap=L_map_val)), 'o-', label='dr')

    ax15.plot(np.diff(colloc_pt), 'o-', color='#83C167', label='dx')
    ax15_1.plot(weights, 'o-', color='m', label=f'LG Quadrature weights (n={L+1})')

    ax17.imshow(S_matrix_real, interpolation=interpolation, cmap='jet')
    ax18.imshow(S_matrix_imag, interpolation=interpolation, cmap='jet')
    ax17_1.imshow(np.abs(S_matrix_real + 1j*S_matrix_imag), interpolation=interpolation, cmap='jet')
    ax17_2.imshow(np.abs(S_matrix_real + 1j*S_matrix_imag)**2, interpolation=interpolation, cmap='jet')


    ax19.plot(r_nm, along_S_diag, 'o-', label=rf'$\sum_k$  ψ$_k$(r) $\cdot$ ψ$_k$(r) ~ l={l}')
    ax19.fill_between(r_nm, along_S_diag, alpha=0.2)

    E_array = [E_field(ti) for ti in t]
    ax20.plot(t/T, E_array, lw=1.5)
    ax20.plot(t[np.argmax(E_array)] / T, max(E_array), 'ro')
    ax20.set_ylim(2*min(E_array), 2*max(E_array))
    ax20.fill_between(t/T, E_array, alpha=0.2)

    max_deformed_pot =  v - E0_au * r
    ax21.plot(r, max_deformed_pot, color='m')
    ax21.fill_between(r, max_deformed_pot, 0, color='m', alpha=0.15)
    for i in range(len(l0_levels)):
        ax21.axhline(l0_levels[i], color='#58C4DD', label='l=0 levels' if i == 0 else '')
    ax21.axis([min(r), 30, -1, 0.1])


    ax1.legend(loc='upper right', ncol=3, columnspacing=0.5, fontsize=12, framealpha=0.8, edgecolor='w')
    ax2.legend(loc='upper right', ncol=3, columnspacing=0.5, fontsize=12, framealpha=0.8, edgecolor='w')
    ax6.legend(loc='upper right', ncol=2, columnspacing=0.5, fontsize=13, framealpha=0.9, edgecolor='w')
    ax3.legend(loc='lower right', fontsize=12, framealpha=0.8, edgecolor='w')
    ax4.legend(loc='lower right', fontsize=12, framealpha=0.8, edgecolor='w')
    ax5.legend(loc='lower right', fontsize=12, framealpha=0.8, edgecolor='w')
    ax7.legend(loc='lower right', fontsize=12, framealpha=0.8, edgecolor='w')
    ax8.legend(loc='upper right', fontsize=15, framealpha=0.8, edgecolor='w')
    ax8a.legend(loc='upper right', fontsize=15, framealpha=0.8, edgecolor='w')
    ax8b.legend(loc='upper right', fontsize=15, framealpha=0.8, edgecolor='w')
    ax9.legend(loc='upper right', fontsize=15, framealpha=0.8, edgecolor='w')
    ax11.legend(loc='lower center', ncol=2, fontsize=12, framealpha=0.8, edgecolor='w')
    ax12.legend(loc='lower right', fontsize=12, framealpha=0.8, edgecolor='w')
    ax13.legend(loc='upper right', fontsize=12, framealpha=0.8, edgecolor='w')
    ax14.legend(loc='upper left', ncol=2, fontsize=12, framealpha=0.8, edgecolor='w')
    ax15.legend(loc='upper left', fontsize=12, framealpha=0.8, edgecolor='w')
    ax15_1.legend(loc='lower center', fontsize=12, framealpha=0.8, edgecolor='w')
    ax16.legend(loc='upper left', fontsize=12, framealpha=0.8, edgecolor='w')
    ax19.legend(loc='best', fontsize=15, framealpha=0.8, edgecolor='w')

    ax7.set_title(f'Potential & Energy levels; l={l}', fontsize=15)
    ax6.set_title(r'Radial functions: |u(r)|$^2$ ; [u(r) = rR(r)]', fontsize=15)

    ax1.set_xlabel('r (a.u) :' + r'$\longrightarrow$', fontsize=15)
    ax2.set_xlabel('r (nm) :' + r'$\longrightarrow$', fontsize=15)
    ax3.set_xlabel('r (nm) :' + r'$\longrightarrow$', fontsize=15)
    ax3.set_ylabel(r'E$_n$ (a.u) :' + r'$\longrightarrow$', fontsize=15)
    ax1.set_ylabel('u(r)' + r'$\longrightarrow$', fontsize=15)
    ax2.set_ylabel(r'|u(r)|$^2$' + r'$\longrightarrow$', fontsize=15)
    ax4.set_xlabel('n :' + r'$\longrightarrow$', fontsize=15)
    ax5.set_xlabel('n :' + r'$\longrightarrow$', fontsize=15)
    ax9.set_xlabel('r (nm) :' + r'$\longrightarrow$', fontsize=15)
    ax6.set_xlabel('r (a.u) :' + r'$\longrightarrow$', fontsize=15)
    ax7.set_xlabel('r (nm) :' + r'$\longrightarrow$', fontsize=15)
    ax7.set_ylabel(r'E$_n$ (a.u) :' + r'$\longrightarrow$', fontsize=15)
    ax8.set_ylabel('u(r) :' + r'$\longrightarrow$', fontsize=15)
    ax8a.set_ylabel('A(x) :' + r'$\longrightarrow$', fontsize=15)
    ax9.set_xlabel('r (nm) :' + r'$\longrightarrow$', fontsize=15)
    ax8b.set_xlabel('x :' + r'$\longrightarrow$', fontsize=15)
    ax9.set_ylabel(r'|u(r)|$^2$ :' + r'$\longrightarrow$', fontsize=15)
    ax8b.set_ylabel(r'|A(x)|$^2$ :' + r'$\longrightarrow$', fontsize=15)
    ax13.set_xlabel('x :' + r'$\longrightarrow$', fontsize=15)
    ax14.set_xlabel(r'collocation points (x$_i$)' + r'$\longrightarrow$', fontsize=15)
    ax14.set_ylabel(r'r$_i$=f(x$_i$)' + r'$\longrightarrow$', fontsize=15)
    ax10.set_xlabel(f'H matrix; (l={l}); interpolated: {interpolation}', fontsize=15)
    ax20.set_xlabel(r't $\longrightarrow$', fontsize=15)
    ax20.set_ylabel(r'E(t) $\longrightarrow$', fontsize=15)
    ax17.set_xlabel(r'$\mathfrak{Re}$(S)', fontsize=20)
    ax18.set_xlabel(r'$\mathfrak{Im}$(S)', fontsize=20)
    ax17_1.set_xlabel(f'|S({l})|', fontsize=20)
    ax17_2.set_xlabel(rf'|S({l})|$^2$', fontsize=20)
    fig6.suptitle(f'S_matrix; (l={l}, k_max={E_n}, dt={dt}); interpolation: {interpolation}', fontsize=20)
    fig6a.suptitle(f'S_matrix; (l={l}, k_max={E_n}, dt={dt}); interpolation: {interpolation}', fontsize=20)


    ax4.set_xticks(ref_n_array)
    ax5.set_xticks(ref_n_array)

    n_values = ["n=" + str(n_ind + l) for n_ind in range(1, len(psi)+1)]
    ax6.set_yticks(offset_arr, n_values, fontsize=12)
    energy_labels = [f'{np.round(value * 27.21, 2)} eV' for value in E]

    ax6_twin.set_yticks(offset_arr, energy_labels, fontsize=12)
    ax6_twin.set_ylim(ax6.get_ylim())
    ax3.set_ylim(ax4.get_ylim())
    ax5.set_ylim(ax4.get_ylim())
    fig1.subplots_adjust(top=0.966, bottom=0.078, left=0.069, right=0.975, wspace=0.196, hspace=0.21)
    fig2.subplots_adjust(top=0.912, bottom=0.097, left=0.061, right=0.934, wspace=0.107)
    fig3.subplots_adjust(top=0.935, bottom=0.09, left=0.109, right=0.928, hspace=0.112)
    fig3a.subplots_adjust(top=0.935, bottom=0.09, left=0.109, right=0.928, hspace=0.112)
    fig4.subplots_adjust(top=0.96, bottom=0.09, left=0, right=0.98, wspace=0.051, hspace=0.16)
    fig5.subplots_adjust(top=0.924, bottom=0.078, left=0.069, right=0.975, wspace=0.116, hspace=0.251)
    fig7.subplots_adjust(top=0.924, bottom=0.078, left=0.069, right=0.975, wspace=0.116, hspace=0.251)
    fig6.subplots_adjust(left=0.02, right=0.98, wspace=0.05)
    fig6a.subplots_adjust(left=0.02, right=0.98, wspace=0.05)

    end_time = time.time()
    print(f'Execution Time (h, m, s) : {secs_to_hr_min_sec(end_time - start_time)}')

    plt.show()
