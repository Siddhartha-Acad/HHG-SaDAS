"""
File: Algo-1_full_analysis.py
Project: HHG-SaDAS
Code Description:
    |---> [Shows detailed analysis of the calculated collocation points]
    |
    | Following Appendix-A of my thesis:
    | - Algo-1 uses a large, equispaced grid in the interval x ∈ (-1, 1)
    | and directly estimates the roots of P'_N(x) from the local maxima
    | of -P'_N(x)² (serving as the initial guesses).

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Collocation points are calculated using the zeros of P'_N (analytical derivative) with fsolve function of scipy.
- It does not save data. Rather, it shows detailed analysis of the calculated values.
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
from scipy.optimize import approx_fprime

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

def P_N(x):
    return legendre(N)(x)

def P_N_deriv(x):
    x = float(x)
    epsilon = dx * 0.01
    return approx_fprime([x], lambda arr: P_N(arr[0]), epsilon)[0]

def P_N_AnaDeriv(x):
    return legendre(N-1)(x) - legendre(N+1)(x)



N = 10
nop_per_segment_values = np.arange(5, 10000, 5, dtype=int)
print(f'~~~~~~~~: Algo-1 :: N = {N} :~~~~~~~~~~')

for nop_per_segment in nop_per_segment_values:
    nopx = nop_per_segment * 10
    total_segments = nopx // nop_per_segment
    x_segment_pts = np.linspace(-1, 1, total_segments + 1)
    dx = float((x_segment_pts[1] - x_segment_pts[0]) / (nop_per_segment - 1))

    print(f'calculating roots for: nop_per_segment={nop_per_segment} (grid_size={nopx})')

    if N % 2 == 0:
        all_roots = [0]
    else:
        all_roots = []

    for i in range(total_segments):
        x_segment = np.linspace(x_segment_pts[i], x_segment_pts[i + 1], nop_per_segment)
        PN_deriv_array = np.array([approx_fprime([float(xi)], lambda arr: P_N(arr[0]), dx)[0] for xi in x_segment])
        pks_at = find_peaks(-PN_deriv_array ** 2)[0]
        colloc_pt_segment = np.array(
            [fsolve(P_N_deriv, x_segment[pks_at[j]], xtol=10 ** -15)[0] for j in range(len(pks_at))]
        )
        all_roots.extend(colloc_pt_segment)

    all_roots = sorted(set(all_roots))
    print(f'    Number of Roots found: {len(all_roots)}')
    if len(all_roots) == N - 1:
        total_grid_size = nop_per_segment * total_segments
        print(f"The minimum grid size required is {total_grid_size}.")
        break



print('\n~~~~~~~~~~: Error Analysis :~~~~~~~~~~')
colloc_pt = all_roots
errors = P_N_AnaDeriv(colloc_pt)
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
print(f'Using: nop_per_segment={nop_per_segment}; grid size={total_grid_size}')

# Measure execution time multiple times (like %%timeit)
num_runs = 7
execution_times = []

for _ in range(num_runs):
    start_time = time.time()

    # Re-run the root-finding with the optimized grid size
    if N % 2 == 0:
        all_roots = [0]
    else:
        all_roots = []

    for i in range(total_segments):
        x_segment = np.linspace(x_segment_pts[i], x_segment_pts[i + 1], nop_per_segment)
        PN_deriv_array = np.array([approx_fprime([float(xi)], lambda arr: P_N(arr[0]), dx)[0] for xi in x_segment])
        pks_at = find_peaks(-PN_deriv_array ** 2)[0]
        colloc_pt_segment = np.array(
            [fsolve(P_N_deriv, x_segment[pks_at[j]], xtol=10 ** -15)[0] for j in range(len(pks_at))]
        )

        all_roots.extend(colloc_pt_segment)

    execution_times.append(time.time() - start_time)

mean_time = np.mean(execution_times)
std_dev_time = np.std(execution_times, ddof=1)  # Sample standard deviation (ddof=1)

# Function to convert time to appropriate unit
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

