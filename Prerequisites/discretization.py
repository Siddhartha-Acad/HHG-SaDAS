"""
File: discretization.py
Project: HHG-SaDAS
Code Description:
    | This script compares two types of radial discretization in polar coordinates.
    | The first subplot shows the old discretization with equispaced radial points.
    | The second subplot shows the new discretization, where the radial coordinate
    | is mapped using a transformation involving parameter L (or L_map). Both
    | discretizations are visualized in 2D using scatter plots.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Old discretization: r is linearly spaced between 0 and 1.
- New discretization: r is transformed as rt = L*(1+x)/(1-x+alpha), giving more
  control over radial density.
- Visualization highlights structural differences between the two discretizations.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import numpy as np
import matplotlib.pyplot as plt
from Assistant.Decorate_axes import decorate_axes_L as da


fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
da.decorate_2d([ax1, ax2], grid=False)

# ~~~~~~~~~~~~~~~~~: Old discretization :~~~~~~~~~~~~~~~~~
r = np.linspace(0, 1, 20)
theta = np.linspace(0, 2*np.pi, 100)

r_m, theta_m = np.meshgrid(r, theta)
x = r_m * np.cos(theta_m)
y = r_m * np.sin(theta_m)
ax1.scatter(x, y, color='indigo')
ax1.set_title('Equispaced radial discretization', fontsize=15)
ax1.axis('square')
# ax1.axis('off')

# ~~~~~~~~~~~~~~~~~: New discretization :~~~~~~~~~~~~~~~~~
x = np.linspace(-1, 1, 50)
L = 0.7
alpha = 2*L / max(r)
rt = L*(1+x) / (1-x+alpha)

rt_m, theta_tm = np.meshgrid(rt, theta)
xt = rt_m * np.cos(theta_tm)
yt = rt_m * np.sin(theta_tm)


print(min(rt), max(rt))
color_values = np.exp(-rt_m**2).flatten()
ax2.scatter(xt, yt, color='indigo')
# ax2.scatter(xt, yt, c=color_values, cmap='jet')           # a bit fun with the colors: as because we have to improvise cmaps to the grid points
ax2.set_title('mapped radial discretization', fontsize=15)
ax2.axis('square')
# ax2.axis('off')
plt.show()

