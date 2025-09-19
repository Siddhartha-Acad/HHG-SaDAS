"""
File: accuracy_colloc_pt.py
Project: HHG-SaDAS
Code Description:
    | This script analyzes the accuracy of collocation points by plugging the collocation points into the analytical form of P'_N(x)
    | It visualizes the results through:
    |   1. A plot of collocation points.
    |   2. A line plot of the analytical error.
    |   3. A bar plot showing the error distribution across all collocation points.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- The script reads collocation points from a text file, calculates the error,
  and generates plots to assess numerical accuracy visually.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.special import legendre
from Assistant.Decorate_axes import decorate_axes_D as da
this_dir = Path(__file__).resolve().parent  # Relative file path system


def a_derivative_P_N(x):
    return legendre(N-1)(x) - legendre(N+1)(x)


N = 14
colloc_file = f'Algo-3_N={N}_AnaDeriv_collocation_points.txt'
colloc_file = this_dir / 'AnaDeriv_Colloc_pt' / colloc_file
colloc_pt = np.loadtxt(colloc_file, skiprows=1, usecols=0)

# colloc_file = 'Algo-3_N=400_AnaDeriv_collocation_points_with_P400.txt'
# colloc_file = this_dir / 'AnaDeriv_Colloc_pt' / colloc_file
# colloc_file = np.loadtxt(colloc_file, skiprows=1)
# colloc_pt = colloc_file[:, 0]
# PN_x   = colloc_file[:, 1]


print('\n~~~~~~~~~~: Error Analysis :~~~~~~~~~~')
errors = a_derivative_P_N(colloc_pt)
mean_error = np.mean(errors)
std_dev_error = np.std(errors, ddof=1)  # Using ddof=1 for sample standard deviation

# for i in range(N-1):
#     print(f'xj[{i}]= {colloc_pt[i]:.5f}; err= {errors[i]:.5e}')

max_abs_value = max(abs(mean_error), abs(std_dev_error))
power = int(np.floor(np.log10(max_abs_value)))
scaled_mean = mean_error / 10**power
scaled_std_dev = std_dev_error / 10**power

print(f'Number of collocation points     : {len(colloc_pt)}')
print(f"Error: mean ± standard deviation : ({scaled_mean:.1f} ± {scaled_std_dev:.1f}) × 10^{power}")


# ---------- Plot ----------
fig, axes = plt.subplots(3, 1, figsize=(10, 12))  # 3 subplots

da.decorate_2d([axes[0]])
axes[0].plot(colloc_pt, 'o-', label='Collocation points')
axes[0].legend(loc='upper left', fontsize=12, framealpha=0.5, edgecolor='k')

# Error as line
da.decorate_2d([axes[1]])
axes[1].plot(errors, 'o-', color='m', label="P'_N(xj)")
axes[1].legend(loc='upper center', fontsize=12, framealpha=0.5, edgecolor='k')

# Error distribution as bar plot
da.decorate_2d([axes[2]])
axes[2].bar(range(len(errors)), errors, color='c', edgecolor='k', alpha=0.7)
axes[2].set_xlabel('Collocation Point Index')
axes[2].set_ylabel('Error')

plt.tight_layout()
plt.show()