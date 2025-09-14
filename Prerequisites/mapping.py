"""
File: mapping.py
Project: HHG-SaDAS
Code Description:
    | This script visualizes the mapping between variable x and its transformed
    | function f(x). The function f(x) is defined with parameters L_map and r_max,
    | and the plot illustrates the mapping as taking points from the horizontal
    | axis (independent variable x) and projecting them to the vertical axis.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- The function f(x) = L_map * (1 + x) / (1 - x + alpha) maps values from x in
  the range [-1, 1] to r.
- The visualization highlights the mapping with dashed vertical and horizontal
  lines, along with markers on the curve.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
from Assistant.Decorate_axes import decorate_axes_L as da

fig = plt.figure()
ax1 = fig.add_subplot(111)
da.decorate_2d(ax1)

def f(x):
    L_map = 0.2
    r_max = 1
    alpha = 2 * L_map / r_max
    return L_map * (1 + x) / (1 - x + alpha)

x = np.linspace(-1, 1, 25)
r = f(x)

ax1.vlines(x, ymin=0, ymax=r, colors='gray', linestyles='dashed', linewidth=1.5)
ax1.hlines(r, xmin=-1, xmax=x, colors='gray', linestyles='dashed', linewidth=1.5)
ax1.plot(x, r, 'o-', c='indigo')


ax1.set_title('Mapping from x to f(x)', fontsize=15)
ax1.set_xlabel(r'x :$\longrightarrow$', fontsize=15)
ax1.set_ylabel(r'r=f(x) :$\longrightarrow$', fontsize=15)
fig.subplots_adjust(top=0.94, bottom=0.069, right=0.98, left=0.048)

plt.show()
