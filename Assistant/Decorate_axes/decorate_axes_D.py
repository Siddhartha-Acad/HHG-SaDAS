"""
File: decorate_axes_D.py
Project: HHG-SaDAS
Code Description: (Dark Theme)
    | This module provides utility functions to format and customize matplotlib
    | axes styles for 2D, 3D, and polar plots. It also supports consistent
    | color schemes, legend formatting, and slider axis styling.

The functions are designed to reduce repetitive code when preparing
publication-quality or presentation-ready figures.

Usage Example
-------------
>>> import matplotlib.pyplot as plt
>>> import Assistant.Decorate_axes.decorate_axes_D as da
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
import matplotlib as mpl
import matplotlib.pyplot as plt
from Assistant.Color_schemes import manim_colors as mc

dec_color = np.concatenate((mc.C, mc.named_color_1, mc.des_col_2))
# dec_color = np.concatenate((mc.C_L, mc.des_col_1))

plt.style.use('dark_background')
plt.rc('font', **{'family': 'serif', 'size': 10})
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=dec_color)


def decorate_imshow(axes_list):
    if not isinstance(axes_list, list):
        axes_list = [axes_list]  # Convert to list if it's not already
    for ax in axes_list:
        ax.grid(False)


def decorate_2d(axes_list, plot_type='line', tick_param=True, grid=True, grid_minorticks=False,
                minor_grid_col='#0078FF', visible_spine='left, bottom', axis_ticks=True, axis='on'):
    # Spine configuration lookup table
    spine_configs = {
        'left, bottom': ['top', 'right'],
        'left, right': ['top'],
        'left': ['top', 'bottom', 'right'],
        'right': ['top', 'bottom', 'left'],
        'top': ['right', 'bottom', 'left'],
        'bottom': ['top', 'right', 'left'],
        'none': ['top', 'bottom', 'left', 'right'],
        'all': []  # Empty list means don't hide any spines
    }

    if not isinstance(axes_list, list):
        axes_list = [axes_list]

    spines_to_hide = spine_configs.get(visible_spine, [])

    for ax in axes_list:
        if not axis_ticks:
            ax.set_xticks([])
            ax.set_yticks([])

        # Hide specified spines
        if spines_to_hide:
            ax.spines[spines_to_hide].set_visible(False)
        elif visible_spine == 'all':
            ax.spines[['top', 'bottom', 'left', 'right']].set_visible(True)

        if grid:
            if grid_minorticks:
                ax.minorticks_on()
                ax.grid(True, alpha=0.5, which='major')
                ax.grid(True, alpha=0.1, color=minor_grid_col, which='minor')
            elif plot_type != 'contourf':
                ax.grid(True, lw=0.4, alpha=0.5, zorder=0)
            else:
                ax.set_aspect('equal', adjustable='box')

        if tick_param:
            ax.tick_params(color='red', width=4)


def da_legend(axes_list, loc='upper right', fontsize=12):
    if not isinstance(axes_list, list):
        axes_list = [axes_list]  # Convert to list if it's not already

    for ax in axes_list:
        ax.legend(loc=loc, fontsize=fontsize, framealpha=0.5, edgecolor='k')


def decorate_polar(axes_list):
    # ax21.set_rlabel_position(45)
    # ax21.set_theta_direction(-1)
    for ax in axes_list:
        ax.grid(True, lw=0.4, alpha=0.5, zorder=0)
        # ax.spines['polar'].set_visible(False)
        ax.axis('off')
        ax.set_aspect('equal')


def decorate_polar_line(axes_list):
    for ax in axes_list:
        ax.grid(True, alpha=0.2, zorder=0)
        ax.spines['polar'].set_visible(False)


def decorate_3d(fig, ax, plot_type='3d_line'):
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.w_xaxis.line.set_color((0, 0, 0, 0))
    ax.w_yaxis.line.set_color((0, 0, 0, 0))
    ax.w_zaxis.line.set_color((0, 0, 0, 0))
    if plot_type != 'surface':
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])


def set_3d_axis(ax, x_min, x_max, y_min, y_max, z_min, z_max):
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)


def decorate_slider(slider_axes):
    slider_axes.spines['top'].set_alpha(0.2)
    slider_axes.spines['right'].set_alpha(0.2)
    slider_axes.spines['left'].set_alpha(0.2)
    slider_axes.spines['bottom'].set_alpha(0.2)

