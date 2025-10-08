import numpy as np
import matplotlib.pyplot as plt
from Assistant.Decorate_axes import decorate_axes_D as da
from scipy.special import jn

data = np.loadtxt('bessel_fwd.dat')

x = 1.0
n_array = data[:, 0]
Jn_rec = data[:, 1]
Jn_f90 = data[:, 2]

rel_err = abs(Jn_f90 - Jn_rec) / Jn_f90

Jn_py = [jn(n, x) for n in n_array]

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
da.decorate_2d([ax1, ax2])

ax1.plot(n_array, Jn_rec, 'o-', markersize=8, label='fwd_recurrence')
ax1.plot(n_array, Jn_f90, 'o-', markersize=6, label='fortran')
ax1.plot(n_array, Jn_py, 'o-', markersize=2, label='python')
ax2.plot(n_array, rel_err, 'o-', label='rel_err')

ax1.set_xticks(n_array)
ax2.set_xticks(n_array)
ax1.legend(fontsize=15, edgecolor='k')
ax2.legend(fontsize=15, edgecolor='k')
plt.tight_layout()
plt.show()