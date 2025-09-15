"""
File: ass_legendre_plot.py
Project: HHG-SaDAS
Code Description:
    | This script visualizes the associated Legendre polynomial P_l^m(cos θ)
    | using `scipy.special.lpmv`. A specific degree `l` and order `m` are chosen,
    | and the function is evaluated over a range of θ values.
    |
    | The plot shows P_l^m(cos θ) against cos(θ), with proper handling of the
    | Condon–Shortley phase factor (-1)^m.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- The `scipy.special.lpmv` function includes the Condon–Shortley phase by default;
  here, it is explicitly corrected with the (-1)^m factor for consistency.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lpmv as a_legendre
import Assistant.Decorate_axes.decorate_axes_D as da

fig = plt.figure()
ax1 = fig.add_subplot(111)
da.decorate_2d(ax1)


l = 8; m = 2

theta = np.linspace(0, 2*np.pi, 200)
P_lm = a_legendre(m, l, np.cos(theta)) * (-1)**m     # Avoiding Condon–Shortley phase (-1)^m factor which is implicit to scipy.special.lpmv


ax1.plot(np.cos(theta), P_lm, 'o-', label='scipy.special.lpmv')
da.da_legend(ax1, fontsize=15)
plt.show()