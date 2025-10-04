"""
File: Gauss_Lobatto.py.py
Project: HHG-SaDAS
Code Description:
    | This script implements and compares different numerical integration methods
    | (Simpson's 1/3 rule, Trapezoidal rule, and Gauss–Lobatto quadrature).
    | The program evaluates the definite integral of test functions, computes
    | relative errors against known analytical results, and prints a comparison
    | table for accuracy analysis.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Uses collocation points from precomputed files to perform Gauss–Lobatto quadrature.
- Integration limits and test functions can be adjusted for different cases.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""


import numpy as np
from pathlib import Path
from scipy.special import erf
from scipy.special import legendre
this_dir = Path(__file__).resolve().parent  # Relative file path system

def integrate(y, x):
    """
    Simpson's 1/3rd method
    :Limit of Integration : Upper and lower limit of independent variable of y.
    :param y: Integrand: x dependent. (odd length)
    :param x: Independent variable on which y variable depends. (odd length)
    :return: Integration of y.
    """
    return (np.abs(x[1] - x[0]) / 3) * (y[0] + 4 * np.sum(y[1:len(y) - 1:2]) + 2 * np.sum(y[2:len(y) - 2:2]) + y[-1])

def integrate_trap(F, x):
    """
    Trapezoidal rule
    :Limit of Integration : Upper and lower limit of independent variable of F.
    :param F: Integrand: x dependent. (any length)
    :param x: Independent variable on which F variable depends. (same length as F)
    :return: Integration of F.
    """
    dx = x[1] - x[0]
    sum = 0
    for i in range(1, len(F) - 1):
        sum += F[i]
    sum = (F[0] + 2 * sum + F[-1]) * dx / 2
    return sum


def f(x): return 1 / (1+x**2)

a, b = 0, 1                              # integration limits
N = 14                                   # number of grid points [Thesis: 4, 8, 14] :-> N+1 grid points

# exact = 2                              # np.sin(x) over [0, np.pi]
# exact = 1 - 0.5*np.sin(2)              # np.sin(x)**2 over [-1, 1]
# exact = 1/28                           # x**27 over [0, 1]
# exact = np.sqrt(np.pi) * erf(1)        # np.exp(-x**2) over [-1, 1]
exact = np.pi / 4               # 1 / (1 + x**2) over [-5, 5]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~: Calculation :~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colloc_file = f'Algo-3_N={N}_Gauss_Lobatto_collocation_points.txt'
colloc_file = this_dir.parent / 'Harmonic_Generation' / 'collocation_points' / 'generator' / colloc_file
colloc_pt = np.loadtxt(colloc_file, skiprows=1, usecols=0)
colloc_pt = np.insert(colloc_pt, [0, len(colloc_pt)], [-1, 1])
int_w = 2 / (N * (N + 1) * (legendre(N)(colloc_pt))**2)

x_mapped = 0.5 * (b - a) * colloc_pt + 0.5 * (b + a)
w_mapped = 0.5 * (b - a) * int_w
integral_GL = sum(w_mapped * f(x_mapped))

x_int = np.linspace(a, b, N+1)
integral_trap = integrate_trap(f(x_int), x_int)
integral_simp = integrate(f(x_int), x_int)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

error_GL = exact - integral_GL
error_simp = exact - integral_simp
error_trap = exact - integral_trap

rel_error_GL = np.abs(error_GL / exact)
rel_error_simp = np.abs(error_simp / exact)
rel_error_trap = np.abs(error_trap / exact)

print('\n~~~ Gauss-Lobatto Quadrature vs Other Methods ~~~')
print(f'Limits         : [{a}, {b}]\n')

print(f'Grid Points    : {N+1}\n')

print(f'Exact Value    : {exact:.15f}')
print(f'Gauss-Lobatto  : {integral_GL:.15f}')
print(f"Simpson's Rule : {integral_simp:.15f}")
print(f'Trapezoidal    : {integral_trap:.15f}\n')

print(f'Rel. Error (Gauss-Lobatto) : {rel_error_GL:.2e}')
print(f"Rel. Error (Simpson's 1/3) : {rel_error_simp:.2e}")
print(f'Rel. Error (Trapezoidal)   : {rel_error_trap:.2e}')

