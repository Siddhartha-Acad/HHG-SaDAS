"""
File: legendre_plot.py
Project: HHG-SaDAS
Code Description:
    | This script visualizes the Legendre polynomial of degree `l` using two
    | representations:
    | 1. Polar-like 2D projection of |P_l(cos θ)| in the x–y plane.
    | 2. Standard line plot of P_l(cos θ) as a function of θ.


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""


import numpy as np
from scipy.special import legendre
from Assistant.Decorate_axes import decorate_axes_D as da
import matplotlib.pyplot as plt


fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(222)
da.decorate_2d([ax1, ax2])

l = 5
theta = np.linspace(0, 2*np.pi, 200)
P_as = legendre(l)(np.cos(theta))

x_plot = abs(P_as) * np.sin(theta)
y_plot = abs(P_as) * np.cos(theta)

ax1.axis('equal')
ax1.plot(x_plot, y_plot, label=r"$P_l^m(cos(\theta))$")
ax2.plot(theta, P_as, label=r"$P_l^m(cos(\theta))$")
ax1.axhline(0, linestyle=' ', label=f'l = {l}')
ax1.set_xlim(-0.8, 0.8)
ax1.legend(fontsize=15, framealpha=0.5, edgecolor='k')
ax2.legend(fontsize=15, framealpha=0.5, edgecolor='k')

fig.subplots_adjust(top=0.925, bottom=0.08, left=0.06, right=0.975, hspace=0, wspace=0.13)

plt.show()