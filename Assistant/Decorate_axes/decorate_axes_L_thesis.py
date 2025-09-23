"""
File: decorate_axes_D.py
Project: HHG-SaDAS
Code Description: (Light Theme used for thesis)
    | This module provides utility functions to format and customize matplotlib
    | axes styles for 2D, 3D, and polar plots. It also supports consistent
    | color schemes, legend formatting, and slider axis styling.

The functions are designed to reduce repetitive code when preparing
publication-quality or presentation-ready figures.

Usage Example
-------------
>>> import matplotlib.pyplot as plt
>>> import Assistant.Decorate_axes.decorate_axes_L_thesis as da
>>> fig, (ax1, ax2) = plt.subplots(2, 1)
>>> da.decorate_2d([ax1, ax2], grid=True, visible_spine='left, bottom')

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
from Assistant.Color_schemes import manim_colors as mc

dec_color = np.concatenate((mc.C_L_reordered, mc.des_col_1))
# dec_color = np.concatenate((mc.C_L, mc.named_color_1, mc.des_col_2))

width = 6.2                         # Width in inches
height = 3                          # Height in inches
fig_scale_factor = 1.5              # big=2 ; medium=1.5; small=1
tickslabel_size = 14
label_fontsize = 15
fig_size = (fig_scale_factor*width, fig_scale_factor*height)

"""
Inkscape writing fontsize = 23
width (in) | height (in) | fig_scale_factor | tickslabel_size (pt) | label_fontsize (pt)| inkscape scale factor | situation
6.2        | 3           | 1.5              | 14                   | 15                 | 120*1                 | Single axes 
6.2        | 3           | 2                | 18                   | 19                 | 120*2, 110*1          | Two axes side by side
6.2        | 0.6         | 2                | 18.5                 | 19                 | 120*1                 | narrow height vertical axes
6.2        | 4           | 2                | 18.5                 | 19                 |                       | 4x2 axes : Wavefunctions
6.2        | 9           | 1.5              | 14                   | 15                 | 120*2                 | full page figure
6.2        | 10          | 1.5              | 14                   | 15                 | 85                    | full page wavefunction
"""


# although tex fontsize is set to 12pt, using 14pt fontsize in figure looks better.
# plt.rc('font', **{'family': 'serif', 'size': tickslabel_size})
# mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=dec_color)

# legend(fontsize=label_fontsize, framealpha=0.7, edgecolor='w')
# ax1.set_xlabel(r'x_label :$\longrightarrow$', fontsize=label_fontsize)
# fig.subplots_adjust(top=0.925, bottom=0.091, left=0.055, right=0.975, hspace=0.245, wspace=0.13)


def decorate_2d(axes_list, plot_type='line', grid=True,
                visible_spine='left, bottom', axis_ticks=True, grid_alpha=0.6):

    # Spine configuration lookup table
    spine_configs = {
        'right, bottom': {'top': False, 'left': False},
        'left, bottom': {'top': False, 'right': False},
        'left, top': {'bottom': False, 'right': False},
        'left, right': {'top': False},
        'left': {'top': False, 'bottom': False, 'right': False},
        'right': {'top': False, 'bottom': False, 'left': False},
        'top': {'right': False, 'bottom': False, 'left': False},
        'bottom': {'top': False, 'right': False, 'left': False},
        'none': {'top': False, 'bottom': False, 'left': False, 'right': False},
        'all': {'top': True, 'bottom': True, 'left': True, 'right': True}
    }

    if not isinstance(axes_list, list):
        axes_list = [axes_list]

    config = spine_configs.get(visible_spine, {})
    for ax in axes_list:
        if not axis_ticks:
            ax.set_xticks([])
            ax.set_yticks([])

        # Apply spine configuration
        for spine, visible in config.items():
            ax.spines[spine].set_visible(visible)

        if grid:
            if plot_type != 'contourf':
                ax.grid(True, lw=0.4, alpha=grid_alpha, zorder=0)
                ax.set_axisbelow(True)
            else:
                ax.set_aspect('equal', adjustable='box')


def add_scale(ax, max_radius, offset=0.05):
    """
    Add a horizontal scale bar beneath the polar plot to indicate radial distances.

    Parameters:
    ax (matplotlib.axes._subplots.AxesSubplot): The polar plot axis.
    max_radius (float): The maximum radius value to be displayed on the scale bar.
    offset (float): Vertical offset to control the distance from the main polar plot.
    """
    # Get the bounding box of the polar axes
    bbox = ax.get_position()

    # Compute dynamic values based on the polar plot's center and rightmost extent
    center_x = bbox.x0 + bbox.width / 2     # Center of the polar plot
    rightmost_x = bbox.x1                   # Rightmost extent

    scale_ax = ax.figure.add_axes([center_x, offset, rightmost_x - center_x, 0])  # Adjusted dynamically
    scale_ax.set_xlim(0, max_radius)
    scale_ax.set_xticks(np.linspace(0, max_radius, num=3))
    scale_ax.set_yticks([])
    scale_ax.spines[['top', 'right', 'left']].set_visible(False)


def decorate_contourf_polar(axes_list, label_padding=10):
    if not isinstance(axes_list, list):
        axes_list = [axes_list]
    for ax in axes_list:
        ax.grid(False); ax.set_rticks([]); ax.set_yticklabels([]); ax.tick_params(pad=label_padding)


def decorate_polar(axes_list):
    # ax21.set_rlabel_position(45)
    # ax21.set_theta_direction(-1)
    for ax in axes_list:
        ax.grid(True, lw=0.4, alpha=0.5, zorder=0)
        # ax.spines['polar'].set_visible(False)
        ax.axis('off')
        ax.set_aspect('equal')