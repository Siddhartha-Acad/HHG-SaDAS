"""
File: absorber_demo.py
Project: HHG-SaDAS
Code Description:
    | Demonstrates the effect of an absorbing layer on a wavefunction attempting
    | to escape a finite simulation box. A right-moving Gaussian wavepacket is
    | used as the test case. The absorbing layer suppresses reflections by
    | gradually damping the wavefunction, preventing unphysical backscattering
    | from the boundaries. For comparison, the free evolution of the wavefunction
    | without the absorber is also shown.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- The normalization of both the initial and time-dependent wavefunction is monitored.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))   # Ensure project root (HHG-SaDAS) is in sys.path

import numpy as np
import matplotlib.pyplot as plt
from Assistant.Color_schemes import manim_colors as mc
from Assistant.Decorate_axes import decorate_axes_D as da
from matplotlib.animation import FuncAnimation

def integrate(y, x):
    """
    Simpson's 1/3rd method
    :Limit of Integration : Upper and lower limit of independent variable of y.
    :param y: Integrand: x dependent.
    :param x: Independent variable on which y variable depends.
    :return: Integration of y.
    """
    return (np.abs(x[1] - x[0]) / 3) * (y[0] + 4 * np.sum(y[1:len(y) - 1:2]) + 2 * np.sum(y[2:len(y) - 2:2]) + y[-1])

def normalize(y, x):
    """
    :param y: Array, to be normalized; x dependent.
    :param x: Independent variable on which y variable depends.
    :return: Normalized y (or), y such that => integral(|y|^2) = 1.
    """
    return np.sqrt(1 / integrate(np.abs(y)**2, x)) * y



N = 50                                                   # Number of terms in Psi(x, t).
nop_x, nop_t = 499, 4001                                 # Number of points in time array; determines the speed of anim.
h_cut, m, a = 1, 0.1, 12                                 # constants.
x, dx = np.linspace(0, a, nop_x, retstep=True)          # x array.
t, dt = np.linspace(0, 1, nop_t, retstep=True)     # time array.


ab_pow = 1/4; Ra = 5
def Absorber_func(x):
    """
    Absorbing layer
    """
    if 0 < x <= Ra:
        return 1
    elif Ra < x < a:
        return np.cos(np.pi * (x - Ra) / (2 * (a - Ra))) ** ab_pow
Absorber = np.array([Absorber_func(x) if Absorber_func(x) is not None else 0 for x in x])

def V(x):
    """
    Infinite potential well
    :param x: x array.
    :return: Potential function.
    """
    inf = 1000
    vx = np.zeros(len(x))
    vx[0], vx[-1] = inf, inf
    return vx

def phi(n, t):
    """
    :param n: Eigenstate number
    :param t: time array.
    :return: Energy of nth energy eigenstate.
    """
    E_n = (n * np.pi * h_cut) ** 2 / (2 * m * a ** 2)
    return np.exp(-1j * E_n * t / h_cut)


# -----------------: INITIAL WAVE FUNCTION :-----------------
def Psi_0(x):
    """
    Initial wave function : GAUSSIAN WAVE FUNCTION
    :param x: x array.
    :return: Initial wave function; Array.
    """
    sigma, mu = 1, 2.3; k = 5
    return (1/np.sqrt(2*np.pi*sigma))*np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) * np.exp(1j*k*x)


fig1 = plt.figure(figsize=(6 * 3.13, 3.1 * 3.13))
ax2 = fig1.add_subplot(111)                           # probability density plot
da.decorate_2d([ax2])

init_Psi = normalize(Psi_0(x), x)                     # Converting the Initial WF into normalized WF
amp = np.max(np.abs(init_Psi) ** 2)                   # Maximum of the normalized WF
ax2.axis([0, a, -0.01, 1.9 * amp])
ax2_inform_cord = 0.02 * a, amp*1.4

def psi_n(n, x):
    """
    :param n: Eigenfunction of nth Energy eigenstate
    :param x: x Array
    :return: Eigenfunction of nth Energy eigenstate
    """
    return np.sqrt(2 / a) * np.sin(n * np.pi * x / a)


def Cn(n):
    """
    :param n: Eigenstate number
    :return: Coefficients of the Superposition of eigenfunctions
    """
    return integrate(psi_n(n, x) * init_Psi, x)


def Psi(x, t):
    """
    NEW PROCESS; FASTER
    :param x: x array
    :param t: time array
    :return: Time evolution function of initial WF
    """
    return sum(Cn(n) * psi_n(n, x) * np.exp(-1j * ((n * np.pi * h_cut) ** 2 / (2 * m * a ** 2)) * t / h_cut) for n in range(1, N))


ax2.plot(x, V(x), color=mc.C[2], lw=1, label='v(x)')
ax2.plot(x, np.abs(init_Psi) ** 2, label='|' + r'$\Psi$' + '(x, 0)|' + r'$^2$', color='red')
ax2.plot(x, Absorber, label='A(x)')
Mod_PSI, = ax2.plot([], [], color='yellow', lw=1.7, label='|' + r'$\Psi$' + '(x, t)|' + r'$^2$')                   # without absorbing function
Mod_PSI_ab, = ax2.plot([], [], color='orange', lw=1.7, label='|A(x)' + r'$\cdot \Psi$' + '(x, t)|' + r'$^2$')      # absorbing layer included

def initialize1():
    Mod_PSI.set_data([], [])
    Mod_PSI_ab.set_data([], [])
    return [Mod_PSI, Mod_PSI_ab, ]

def update1(i):
    PSI = Psi(x, t[i])
    mod_PSI_sq = (np.abs(PSI)) ** 2
    mod_PSI_sq_ab = (np.abs(Absorber * PSI)) ** 2

    Mod_PSI.set_data(x, mod_PSI_sq)
    Mod_PSI_ab.set_data(x, mod_PSI_sq_ab)

    str_time = 't = ' + str("{:.5f}".format(t[i])) + '; step = ' + str(i) + '\n' + r'$\int_{0}^{a}$' + '|' + r'$\Psi$' + \
               '(x, t)' + '|' + r'$^2$' + 'dx' + ' = ' + str(np.round(integrate(mod_PSI_sq, x), 5)) +\
               '\n' + r'$\int_{0}^{a}$' + '|' + r'$\Psi$' + '(x, 0)' + '|' + r'$^2$' + 'dx' + ' = ' + \
               str(np.round(integrate(np.abs(init_Psi) ** 2, x), 5))
    annotate_time = ax2.annotate(str_time, xy=ax2_inform_cord, fontsize=15)
    return [Mod_PSI, Mod_PSI_ab, annotate_time, ]

anim1 = FuncAnimation(fig1, update1, frames=len(t), init_func=initialize1, interval=1, blit=True)


ax2.set_title('|' + r'$\Psi$' + '(x, t)|' + r'$^2$' + f' ; N = {N}', fontsize=20)
ax2.set_xlabel('X :' + r'$\longrightarrow$', fontsize=15)
ax2.legend(loc='upper right', fontsize=15, frameon=False)
fig1.suptitle('Gaussian wavepacket, bouncing between walls', color=mc.C[2], fontsize=22)
fig1.subplots_adjust(bottom=0.083, top=0.886, left=0.057, right=0.962, wspace=0.126)

plt.show()