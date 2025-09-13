"""
File: 1d_func_integ.py
Project: HHG-SaDAS
Code Description:
    | This script demonstrates numerical integration of f(x) over [a,b] using
    | three different methods:
    | 1. Gauss–Legendre Quadrature
    | 2. Simpson’s 1/3rd Rule
    | 3. Trapezoidal Rule
    |
    | It compares the computed values with the exact integral and reports the
    | relative error for Gauss–Legendre Quadrature.


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Simpson’s and Trapezoidal methods are implemented manually for comparison.
- This code provides a simple validation of quadrature rules against known
  analytical results.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas
  Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import numpy as np

def f(x):
    return np.sin(x)

def integrate_trap(F, x):
    dx = x[1] - x[0]
    sum = 0
    for i in range(1, len(F) - 1):
        sum += F[i]
    return (F[0] + 2 * sum + F[-1]) * dx / 2

def integrate(y, x):
    """
    Simpson's 1/3rd method
    :Limit of Integration : Upper and lower limit of independent variable of y.
    :param y: Integrand: x dependent.
    :param x: Independent variable on which y variable depends.
    :return: Integration of y.
    """
    return (np.abs(x[1] - x[0]) / 3) * (y[0] + 4 * np.sum(y[1:len(y) - 1:2]) + 2 * np.sum(y[2:len(y) - 2:2]) + y[-1])

n = 5
a, b = 0, 1
x, w = np.polynomial.legendre.leggauss(n)
x_mapped = 0.5 * (b - a) * x + 0.5 * (b + a)
w_mapped = 0.5 * (b - a) * w
integral_GL = sum(w_mapped * f(x_mapped))


# comparison with other method
x_int = np.linspace(a, b, n)
integral_trap = integrate_trap(f(x_int), x_int)
integral_simp = integrate(f(x_int), x_int)          # number of grid points for simpson's 1/3 must be odd.


exact = 1 - np.cos(1)
rel_error = abs(exact - integral_GL) / exact
print('~~~: Gauss Legendre Quadrature :~~~')
print('Integral      : sin(x)')
print(f'limits        : ({a}, {b})\n')
print('nop           :', n)
print('exact value   :', exact)
print('Gauss_Legendre:', integral_GL)
print("Simpson's     :", integral_simp)
print('Trapezoidal   :', integral_trap)
print('Relative Error:', rel_error)