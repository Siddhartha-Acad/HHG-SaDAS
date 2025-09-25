"""
File: Algo-3_full_analysis.py
Project: HHG-SaDAS
Code Description:
    |---> [Shows detailed analysis of the calculated collocation points]
    |
    | Following Appendix-A of my thesis:
    | - Algo-3 uses the non-equispaced grid x⁺(ξ) ∈ (0, 1), as in Eq.A.11,
    |   to determine the roots of Λ_N(x) from the local maxima of -Λ_N(x)² (serving as the initial guesses).
    | - It calculates only half of the roots (those in the positive interval), while the other half
    |   (in the negative interval) are obtained using the parity relation in Eq.A.10.
    | - This algorithm is graphically presented in the flowchart of Fig.A.4.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Collocation points are calculated using the zeros of P'_N (analytical derivative) with fsolve function of scipy.
- Optionally writes collocation points or collocation points with P_N values to a text file.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""
import time
import warnings
import numpy as np
from scipy.optimize import fsolve
from scipy.special import legendre
from scipy.signal import find_peaks
warnings.filterwarnings("ignore", category=RuntimeWarning)

def P_N_AnaDeriv(x, N):
    return legendre(N-1)(x) - legendre(N+1)(x)

def f_rev(x_array):
    r_max = 1
    L_map = 0.2
    alpha = 2 * L_map / r_max
    map_func = L_map * (1 + x_array) / (1 - x_array + alpha)
    return -map_func[::-1] + r_max


N = 200
nopx_values = np.arange(5, 10000, 5, dtype=int)
print(f'~~~~~~~~: Algo-3 :: N = {N} :~~~~~~~~~~')

for nopx in nopx_values:
    print(f'calculating roots for: grid_size={nopx}')

    nopx_half = int(nopx / 2)
    xi = np.linspace(-1, 1, nopx_half)
    x_mapped = f_rev(xi)  # Apply the grid mapping

    PN_deriv_array = P_N_AnaDeriv(x_mapped, N)
    pks_at = find_peaks(-PN_deriv_array ** 2)[0]
    colloc_pt = np.array([fsolve(P_N_AnaDeriv, x_mapped[pks_at[j]], args=(N,))[0] for j in range(len(pks_at))])

    if N % 2 == 0:    # Symmetrizing the roots
        colloc_pt = np.concatenate((-colloc_pt[::-1], [0.0], colloc_pt))
    else:
        colloc_pt = np.concatenate((-colloc_pt[::-1], colloc_pt))
    print(f'    Number of roots found: {len(colloc_pt)}')

    if len(colloc_pt) == N - 1:
        print(f"Minimum grid size required: nopx = {nopx}")
        break


print('\n~~~~~~~~~~: Error Analysis :~~~~~~~~~~')
errors = P_N_AnaDeriv(colloc_pt, N)
mean_error = np.mean(errors)
std_dev_error = np.std(errors, ddof=1)  # Using ddof=1 for sample standard deviation

for i in range(N-1):
    print(f'xj[{i}]= {colloc_pt[i]:.5f}; err= {errors[i]:.5e}')

max_abs_value = max(abs(mean_error), abs(std_dev_error))
power = int(np.floor(np.log10(max_abs_value)))
scaled_mean = mean_error / 10**power
scaled_std_dev = std_dev_error / 10**power

print(f'Number of collocation points     : {len(colloc_pt)}')
print(f"Error: mean ± standard deviation : ({scaled_mean:.1f} ± {scaled_std_dev:.1f}) × 10^{power}")
# print(rf'copy: \( ({scaled_mean:.1f} \pm {scaled_std_dev:.1f})\times 10^{ {power} } \)')


print('\n~~~~~: Measuring Execution time :~~~~~')
print(f'Using: grid size={nopx}')
num_runs = 7; execution_times = []
for _ in range(num_runs):
    start_time = time.time()

    nopx_half = int(nopx / 2)
    xi = np.linspace(-1, 1, nopx_half)
    x_mapped = f_rev(xi)

    PN_deriv_array = P_N_AnaDeriv(x_mapped, N)
    pks_at = find_peaks(-PN_deriv_array ** 2)[0]
    colloc_pt = np.array([fsolve(P_N_AnaDeriv, x_mapped[pks_at[j]], args=(N,))[0] for j in range(len(pks_at))])

    # Symmetrizing the roots
    if N % 2 == 0:
        colloc_pt = np.concatenate((-colloc_pt[::-1], [0.0], colloc_pt))
    else:
        colloc_pt = np.concatenate((-colloc_pt[::-1], colloc_pt))

    execution_times.append(time.time() - start_time)

mean_time = np.mean(execution_times)
std_dev_time = np.std(execution_times, ddof=1)  # Sample standard deviation (ddof=1)

def format_time(value):
    if value < 1e-3:
        return f"{value * 1e6:.1f} µs"
    elif value < 1:
        return f"{value * 1e3:.1f} ms"
    elif value < 60:
        return f"{value:.1f} s"
    elif value < 3600:
        return f"{value / 60:.1f} min"
    else:
        return f"{value / 3600:.1f} hr"

print(f"Execution time: {format_time(mean_time)} ± {format_time(std_dev_time)} per loop "
      f"(mean ± std. dev. of {num_runs} runs, 1 loop each)")
# print(f'copy: \(( {format_time(mean_time)} \pm {format_time(std_dev_time)} )\) ms\n')






