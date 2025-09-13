"""
File: Truncate_cmap.py
Project: HHG-SaDAS

Code Description:
    | This script defines a utility function to truncate an existing
    | Matplotlib colormap to a specified subrange. This is useful when
    | only a portion of a colormap is desired for highlighting specific
    | data ranges.
    |
    | Features:
    |   • Accepts any Matplotlib colormap by name
    |   • Defines the minimum and maximum normalized values to keep
    |   • Returns a new LinearSegmentedColormap of the truncated range
    |   • Supports specifying the number of discrete points in the new colormap


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
-
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""


import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

def truncate_cmap(cmap_in='jet', minval=0.0, maxval=1.0, n=100):
    cmap_in = plt.get_cmap(cmap_in)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_in.name, a=minval, b=maxval),
        cmap_in(np.linspace(minval, maxval, n)))

    return new_cmap
