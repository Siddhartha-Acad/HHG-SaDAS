"""
File: Working_grid_plot.py
Project: HHG-SaDAS
Code Description:
    | Demonstrates a 2D grid with boundary conditions. The interior values are set
    | to one, and the outer boundary is set to zero, then visualized with a color mesh.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- This is only for visualization.
- The actual grid is not rectangular but defined in spherical polar coordinates.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt

N = 20
matrix = np.ones((N, N))
matrix[0, :] = 0
matrix[-1, :] = 0
matrix[:, 0] = 0
matrix[:, -1] = 0

fig, ax = plt.subplots(figsize=(6, 6))
ax.pcolormesh(matrix, cmap='ocean', edgecolors='k')
ax.set_aspect('equal')
plt.show()
