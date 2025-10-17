import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
from Assistant.Decorate_axes import decorate_axes_D as da

data = np.loadtxt('legendre_poly.dat')

n_array = data[:, 0]
Pn_f90 = data[:, 1]

x = 1.0
Pn_scipy = np.array([legendre(n_val)(x) for n_val in n_array])

rel_err = abs(Pn_scipy - Pn_f90) / Pn_scipy

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
da.decorate_2d([ax1, ax2])

ax1.plot(n_array, Pn_f90, 'o-', markersize=8, label='fortran')
ax1.plot(n_array, Pn_scipy, 'o-', markersize=4, label='python')
ax2.plot(n_array, rel_err, 'o-', label='rel_err')

ax1.legend(fontsize=15, edgecolor='k')
plt.tight_layout()
plt.show()