"""
File: cmap_test.py.py
Project: HHG-SaDAS

Code Description:
    | This code demonstrates how to use a custom colormap previously saved
    | as a NumPy array ('my_nipy_spectral.npy') for plotting in Matplotlib.
    | The script:
    |   • Loads the custom colormap data
    |   • Creates a LinearSegmentedColormap from it
    |   • Registers it for direct use in Matplotlib
    |   • Uses it to generate a contour plot of a 2D potential function


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
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
Z = -1 / np.sqrt(X**2 + Y**2)

# Load the saved colormap data
loaded_cmap_data = np.load("my_nipy_spectral.npy")

# Create a new colormap from loaded data
my_nipy_spectral = mcolors.LinearSegmentedColormap.from_list("my_nipy_spectral", loaded_cmap_data)

# Register it in Matplotlib for direct use
plt.register_cmap(cmap=my_nipy_spectral)

# Now you can use it like a built-in colormap!
plt.figure()
plt.contourf(X, Y, Z, levels=500, cmap="my_nipy_spectral")
plt.axis('equal')
plt.axis('square')
plt.colorbar()
plt.show()
