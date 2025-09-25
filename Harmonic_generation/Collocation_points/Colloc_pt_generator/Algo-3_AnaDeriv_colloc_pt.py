"""
File: Algo-3_AnaDeriv_colloc_pt.py
Project: HHG-SaDAS
Code Description:
    | *** [Main Gauss-Lobatto collocation point generating code] ***
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
import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import warnings
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import legendre
from scipy.signal import find_peaks
from Assistant.Decorate_axes import decorate_axes_D as da

this_dir = Path(__file__).resolve().parent                                              # Relative file path system
warnings.filterwarnings("ignore", category=RuntimeWarning)


def f_rev(x_array):             # f(x) reversed to have dense grid towards far.
    r_max = 1; L_map = 0.5
    alpha = 2 * L_map / r_max
    map_func = L_map*(1 + x_array) / (1 - x_array + alpha)
    return -map_func[::-1] + r_max

def P_N_AnaDeriv(x, N):
    return legendre(N-1)(x) - legendre(N+1)(x)


N = 200
Write_PN = False


nopx = 2760
xi = np.linspace(-1, 1, nopx); x_mapped = f_rev(xi)
PN_deriv_array = P_N_AnaDeriv(x_mapped, N)
pks_at = find_peaks(-PN_deriv_array ** 2)[0]
colloc_pt = np.array([fsolve(P_N_AnaDeriv, x_mapped[pks_at[j]], args=(N,))[0] for j in range(len(pks_at))])

if N % 2 == 0:
    colloc_pt = np.concatenate((-colloc_pt[::-1], [0.0], colloc_pt))
else:
    colloc_pt = np.concatenate((-colloc_pt[::-1], colloc_pt))



errors = P_N_AnaDeriv(colloc_pt, N)
mean_error = np.mean(errors)
std_dev_error = np.std(errors, ddof=1)  # Using ddof=1 for sample standard deviation

# for i in range(N-1):
#     print(f'xj[{i}]= {colloc_pt[i]:.5f}; err= {errors[i]:.5e}')

max_abs_value = max(abs(mean_error), abs(std_dev_error))
power = int(np.floor(np.log10(max_abs_value)))

scaled_mean = mean_error / 10 ** power
scaled_std_dev = std_dev_error / 10 ** power
print(f'~~~~~~~~: Algo-3 :: N = {N} :~~~~~~~~~~')
print(f'no. of collocation point  : {len(colloc_pt)}')
print(f"mean ± standard deviation : ({scaled_mean:.2f} ± {scaled_std_dev:.2f}) × 10^{power}")

# Text file introduced on 2/12/2024: listening Pandit Ravi Shankar ESCONDITO 2011 full concert
if len(colloc_pt) == N - 1:
    if not Write_PN:
        txt_file_name = f'Algo-3_{N=}_AnaDeriv_collocation_points.txt'
        file_path = this_dir / txt_file_name
        # if file_path.exists():
        #     raise FileExistsError(f"File already exists: {file_path.name}")

        header = f'{"colloc_pt":>12} (N={N})'
        np.savetxt(file_path, colloc_pt, fmt='%20.15f', header=header, comments='')

    else:
        txt_file_name = f'Algo-3_{N=}_AnaDeriv_collocation_points_with_P{N}.txt'
        file_path = this_dir / txt_file_name
        if file_path.exists():
            raise FileExistsError(f"File already exists: {file_path.name}")

        colloc_pt_and_PN = np.column_stack((colloc_pt, legendre(N)(colloc_pt)))
        header = f'{"colloc_pt":>20} {"P_" + str(N) + "(colloc_pt)":>20}'
        np.savetxt(file_path, colloc_pt_and_PN, fmt='%20.15f', header=header, comments='')

    print(f'file created: {txt_file_name}')


fig = plt.figure()
ax1 = fig.add_subplot(111)
da.decorate_2d(ax1)

ax1.plot(x_mapped, PN_deriv_array, 'o-', markersize=2, color=da.mc.sweet_green, zorder=1)
ax1.scatter(x_mapped[pks_at], np.zeros_like(pks_at), color='w', zorder=2)
plt.show()
