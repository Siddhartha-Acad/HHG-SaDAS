"""
File: Sigmoid_cmap.py
Project: HHG-SaDAS

Code Description:
    | This script defines a function to create a sigmoid-transformed version
    | of any Matplotlib colormap. The sigmoid transformation allows for
    | emphasizing a particular range of values in the colormap.
    |
    | Features:
    |   • Defines a sigmoid function for smooth transition
    |   • Applies the sigmoid to the input colormap
    |   • Returns a new LinearSegmentedColormap
    |   • Includes an example comparison between the original and sigmoid
    |     colormaps using a 10x10 random dataset


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
from matplotlib.colors import LinearSegmentedColormap

x0 = 0.5; a = 0.05

def sigmoid(x, x0, a):
    return 1 / (1 + np.exp(-(x - x0) / a))


def sigmoid_cmap(color_map='Blues', central_position=x0, transition_width=a, N=256):
    original_cmap = getattr(plt.cm, color_map)
    sigmoid_values = sigmoid(np.linspace(0, 1, N), central_position, transition_width)
    sigmoid_colors = original_cmap(sigmoid_values)
    return LinearSegmentedColormap.from_list("SigmoidBlues", sigmoid_colors, N=N)



if __name__ == '__main__':
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    data = np.random.rand(10, 10)
    plot1 = ax1.imshow(data, cmap='Greens')
    fig.colorbar(plot1, ax=ax1, fraction=0.046, pad=0.04)
    ax1.set_title('Normal Cmap', fontsize=15)

    plot2 = ax2.imshow(data, cmap=sigmoid_cmap(color_map='Greens', central_position=0.6, transition_width=0.01))
    fig.colorbar(plot2, ax=ax2, fraction=0.046, pad=0.04)
    ax2.set_title('Sigmoid Cmap', fontsize=15)

    plt.show()