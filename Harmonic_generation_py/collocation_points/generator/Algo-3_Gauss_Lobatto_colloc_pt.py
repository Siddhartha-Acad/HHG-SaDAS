"""
File: Algo-3_Gauss_Lobatto_colloc_pt.py
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

import time
import warnings
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import legendre
from scipy.signal import find_peaks
from Assistant.Decorate_axes import decorate_axes_L as da

this_dir = Path(__file__).resolve().parent                       # Relative file path system
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ~~~~~~~~~~~~~~~~~~~~~~: Common Figure Settings :~~~~~~~~~~~~~~~~~~~~
width = 6.2                         # Width in inches
height = 4                          # Height in inches
fig_scale_factor = 1.5              # big=2 ; medium=1.5; small=1
fig_size = (fig_scale_factor*width, fig_scale_factor*height)

plt.rc('font', **{'family': 'serif', 'size': 14})
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Useful Functions                          |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def f_rev(xi):             # f(x) reversed to have dense grid towards far.
    """
    Reversed nonlinear mapping function.

    Creates a grid denser toward x -> +1
        f(xi) = 1 - L_map * (1 - xi) / (1 + xi + α)
    with α = 2 * L_map.

    Reference
    ----------
    Section A.2.3 — Nonlinear grid discretization and the root-finding algorithm
    (See Equation A.11 and Figure A.3)
    """
    L_map = 0.5
    alpha = 2 * L_map
    return 1.0 - L_map * (1 - xi) / (1 + xi + alpha)


def Lambda(x, N):
    """
    The capital Lambda function Λ_N(x)

    Defined as: Λ_N(x) = P_{N−1}(x) − P_{N+1}(x)

    where P_k(x) is the Legendre polynomial of degree k.

    Reference
    -----------
    Appendix A: "An efficient algorithm to calculate the Gauss–Lobatto collocation points"
    (See Equation A.6 and Figure A.2)
    """
    return legendre(N-1)(x) - legendre(N+1)(x)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Main control                            |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = 200                     # Collocation points (xj) = P'_N(xj) = 0 :: total collocation points = N-1. (my thesis work was done with N=200)
Write_PN = False            # setting True : creates data file (.txt) with P_N(xj) to the second column.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          MAIN computation                          |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start_time = time.perf_counter()

xi = np.linspace(-1, 1, 100000); xi_mapped = f_rev(xi)
PN_deriv_array = Lambda(xi_mapped, N)
pks_at = find_peaks(-PN_deriv_array ** 2)[0]            # The initial guesses
colloc_pt = np.array([                                  # calculating using Newton-Raphson method.
    fsolve(Lambda, xi_mapped[p], args=(N,))[0]
    for p in pks_at
])

if N % 2 == 0:    # using parity to mirror and make complete set of collocation points.
    colloc_pt = np.concatenate((-colloc_pt[::-1], [0.0], colloc_pt))
else:
    colloc_pt = np.concatenate((-colloc_pt[::-1], colloc_pt))

end_time = time.perf_counter()
CPU_time = end_time - start_time


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           Error Analysis                           |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
errors = Lambda(colloc_pt, N)
mean_error = np.mean(errors)
std_dev_error = np.std(errors, ddof=1)  # Using ddof=1 for sample standard deviation

# for i in range(N-1):
#     print(f'xj[{i}]= {colloc_pt[i]:.5f}; err= {errors[i]:.5e}')

max_abs_value = max(abs(mean_error), abs(std_dev_error))
power = int(np.floor(np.log10(max_abs_value)))

scaled_mean = mean_error / 10 ** power
scaled_std_dev = std_dev_error / 10 ** power
print(f'\n~~~~~~~~~~~~~~: Algo-3 :: N = {N} :~~~~~~~~~~~~~~')
print(f'no. of collocation point  : {len(colloc_pt)}')
print(f"mean ± standard deviation : ({scaled_mean:.2f} ± {scaled_std_dev:.2f}) × 10^{power}\n")
print(f'CPU time : {CPU_time:.4f} seconds')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Plotting                              |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig = plt.figure(figsize=fig_size)
ax1 = fig.add_subplot(111)
da.decorate_2d(ax1)

ax1.plot(xi_mapped, PN_deriv_array, lw=2.0, color=da.mc.sweet_green, label=r"P'$_N$(x)", zorder=1)
ax1.scatter(xi_mapped[pks_at], np.zeros_like(pks_at), color='darkviolet', label=r"P'$_N$(x$_j$)", zorder=2)

ax1.set_xlabel(r'$x = f_{rev}(\xi); \xi \in (-1, 1)$', fontsize=15)
ax1.set_ylabel(r"P'$_N$(x)", fontsize=15)
ax1.set_title(f"Algorithm 3: Gauss–Lobatto collocation points ($x_j$ : $P'_N(x_j) = 0$) ; N={N}", fontsize=14, pad=15)
ax1.legend(loc='upper right', fontsize=14, framealpha=0.5, edgecolor='w')
fig.subplots_adjust(
    top=0.905,
    bottom=0.12,
    left=0.125,
    right=0.95,
    hspace=0.2,
    wspace=0.2
)
plt.show()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     Writing data to .txt file                      |
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Text file introduced on 2/12/2024: listening Pandit Ravi Shankar ESCONDITO 2011 full concert
if len(colloc_pt) == N - 1:
    if not Write_PN:
        file_name = f'Algo-3_{N=}_Gauss_Lobatto_collocation_points.txt'
        file_path = this_dir / file_name
        if file_path.exists():                  # Prevents from overwriting existing data file
            raise FileExistsError(f"File already exists: {file_path.name}")

        header = f'{"colloc_pt":>12} (N={N})'
        np.savetxt(file_path, colloc_pt, fmt='%20.15f', header=header, comments='')

    else:
        file_name = f'Algo-3_{N=}_Gauss_Lobatto_collocation_points_with_P{N}.txt'
        file_path = this_dir / file_name
        if file_path.exists():                  # Prevents from overwriting existing data file
            raise FileExistsError(f"File already exists: {file_path.name}")

        colloc_pt_and_PN = np.column_stack((colloc_pt, legendre(N)(colloc_pt)))
        header = f'{"colloc_pt":>20} {"P_" + str(N) + "(colloc_pt)":>20}'
        np.savetxt(file_path, colloc_pt_and_PN, fmt='%20.15f', header=header, comments='')

    print(f'file created: {file_name}')
