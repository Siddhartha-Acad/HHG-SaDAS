"""
File: Algo-2_full_analysis.py
Project: HHG-SaDAS
Code Description:
    |---> [Shows detailed analysis of the calculated collocation points]
    |
    | Following Appendix-A of my thesis:
    | - Algo-2 also uses an equispaced grid in x ∈ (-1, 1). However, unlike Algo-1,
    | it estimates the roots of Λ_N(x), requiring a much smaller grid size to generate
    | the initial guesses.

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
import numpy as np
from scipy.optimize import fsolve
from scipy.special import legendre
from scipy.signal import find_peaks

def P_N_AnaDeriv(x, N):
    return legendre(N-1)(x) - legendre(N+1)(x)

N = 10
nopx_values = np.arange(5, 10000, 5, dtype=int)
print(f'~~~~~~~~: Algo-2 :: N = {N} :~~~~~~~~~~')

for nopx in nopx_values:
    print(f'calculating roots for: grid_size={nopx}')

    x = np.linspace(-1, 1, nopx)
    PN_deriv_array = P_N_AnaDeriv(x, N)
    pks_at = find_peaks(-PN_deriv_array ** 2)[0]
    colloc_pt = np.array([fsolve(P_N_AnaDeriv, x[pks_at[j]], args=(N,))[0] for j in range(len(pks_at))])

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

    x = np.linspace(-1, 1, nopx)
    PN_deriv_array = P_N_AnaDeriv(x, N)
    pks_at = find_peaks(-PN_deriv_array ** 2)[0]
    colloc_pt = np.array([fsolve(P_N_AnaDeriv, x[pks_at[j]], args=(N,))[0] for j in range(len(pks_at))])

    execution_times.append(time.time() - start_time)

mean_time = np.mean(execution_times)
std_dev_time = np.std(execution_times, ddof=1)  # Sample standard deviation (ddof=1)

# Function to format time in appropriate units
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
