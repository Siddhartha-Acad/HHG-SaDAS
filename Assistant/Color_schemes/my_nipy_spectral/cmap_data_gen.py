"""
File: cmap_data_gen.py
Project: HHG-SaDAS

Code Description:
    | This code creates a modified version of the 'nipy_spectral' colormap
    | from Matplotlib, where the highest color values smoothly fade to white.
    | The resulting colormap is saved as a NumPy array for later use in plots
    | within the HHG-SaDAS project.


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

# Get the original colormap
cmap = plt.cm.nipy_spectral
new_cmap = cmap(np.linspace(0, 1, 256))  # Get colors

# Define transition range for smooth fade to white
gray_start_idx = 230  # Where gray starts
white_end_idx = 255  # Where white fully takes over

for i in range(gray_start_idx, white_end_idx + 1):
    blend_factor = (i - gray_start_idx) / (white_end_idx - gray_start_idx)
    new_cmap[i] = new_cmap[i] * (1 - blend_factor) + np.array([1, 1, 1, 1]) * blend_factor

# Save colormap data
np.save("my_nipy_spectral.npy", new_cmap)
print("Colormap saved: my_nipy_spectral.npy")
