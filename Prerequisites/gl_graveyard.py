"""
File: gl_computing_speed.ipynb
Project: HHG-SaDAS
Code Description:
    | Dumpyard of partial wave calculating functions.
    | This code does nothing functionally, but serves to visualize how the partial waves looked
    | when they were developed. It contains different implementations of computing gl(r) for a
    | given wavefunction Psi_t, emphasizing both clarity and computational efficiency.
    |
    | Better shown in: HHG-SaDAS\Prerequisites\gl_museum

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

def psi_t_call(Psi_t, i, j):
    """
    returns ψ(r_i, θ_j, t)
    """
    return Psi_t[j][i]


def gl_1(Psi_t):
    """
    Slower gl(r) calculating
    Calculate partial waves for a given wavefunction.
    Psi_t : Complex array representing the wavefunction at each combination of (r, θ, t).
    Notes : θ must include values θ_k, as they are required for the calculation.
    """
    gl_arr = []
    for l_ind in range(L):                  # this L should come from H_ij_S_ij.py
        g_l = np.zeros(len(r), dtype=np.complex128)
        for i in range(len(g_l)):
            g_k = 0
            for k in range(len(weights)):
                g_k += weights[k] * legendre(l_ind)(roots[k]) * psi_t_call(Psi_t, i, theta_k_locs[k])
            g_l[i] = g_k
        gl_arr.append(g_l * (l_ind+0.5))
    return gl_arr


def gl_2(Psi_t):
    """
    Faster gl(r) calculating
    Calculate partial waves for a given wavefunction.
    Psi_t : Complex array representing the wavefunction at each combination of (r, θ, t).
    Notes : θ must include values θ_k, as they are required for the calculation.
    """
    gl_arr = []
    len_r = len(r); len_k = len(weights)
    legendre_vals = np.array([[legendre(l_index)(root) for root in roots] for l_index in range(L)])
    for l_ind in range(L):
        g_l = np.zeros(len_r, dtype=np.complex128)
        psi_vals = np.array([[psi_t_call(Psi_t, i, theta_k_locs[k]) for k in range(len_k)] for i in range(len_r)])
        for i in range(len_r):
            g_l[i] = np.sum(weights * legendre_vals[l_ind, :] * psi_vals[i, :])
        gl_arr.append(g_l * (l_ind + 0.5))
    return gl_arr


# Precompute Legendre polynomial values
legendre_vals_ext = np.array([[legendre(l_index)(root) for root in roots] for l_index in range(L)])
def gl_3(Psi_t):
    """
    Calculate partial waves for a given wavefunction.
    Psi_t : Complex array representing the wavefunction at each combination of (r, θ, t).
    Notes : θ must include θ_k, as they are required for the calculation.
    ~~ Psi_t[theta_k_locs[k]][i] = ψ(r_i, θ_k, t)
    """
    gl_arr = []; len_r = len(r); len_k = len(weights)
    psi_vals = np.array([[Psi_t[theta_k_locs[k]][i] for k in range(len_k)] for i in range(len_r)])
    for l_ind in range(L):
        g_l = (weights * legendre_vals_ext[l_ind, :]) @ psi_vals.T
        gl_arr.append(g_l * (l_ind + 0.5))
    return gl_arr