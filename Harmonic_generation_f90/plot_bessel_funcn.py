import numpy as np
import matplotlib.pyplot as plt
from Assistant.Decorate_axes import decorate_axes_D as da
from scipy.special import jn

data = np.loadtxt('bessel_funcn.dat')

x = 5.0
n_array = data[:, 0]
Jn_fwd = data[:, 1]
Jn_bwd = data[:, 2]
Jn_f90 = data[:, 3]
Jn_py = [jn(n, x) for n in n_array]


rel_err_fwd = abs(Jn_f90 - Jn_fwd) / Jn_f90
rel_err_bwd = abs(Jn_f90 - Jn_bwd) / Jn_f90


fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
da.decorate_2d([ax1, ax2, ax3, ax4])

ax1.plot(n_array, Jn_fwd, 'o-', markersize=8, label='forward')
ax1.plot(n_array, Jn_f90, 'o-', markersize=6, label='fortran')
ax1.plot(n_array, Jn_py, 'o-', markersize=2, label='python')

ax2.plot(n_array, Jn_bwd, 'o-', markersize=8, label='backward')
ax2.plot(n_array, Jn_f90, 'o-', markersize=6, label='fortran')
ax2.plot(n_array, Jn_py, 'o-', markersize=2, label='python')

ax3.plot(n_array, rel_err_fwd, 'o-', label='rel_err_fwd')
ax4.plot(n_array, rel_err_bwd, 'o-', label='rel_err_bwd')


ax1.set_xlabel('N values', fontsize=15); ax2.set_xlabel('N values', fontsize=15)
ax3.set_xlabel('N values', fontsize=15); ax4.set_xlabel('N values', fontsize=15)
ax1.set_ylabel(rf'J$_{{N}}(x={x:.2f})$', fontsize=15); ax2.set_ylabel(rf'J$_{{N}}(x={x:.2f})$', fontsize=15)
ax3.set_ylabel(rf'J$_{{N}}(x={x:.2f})$', fontsize=15); ax4.set_ylabel(rf'J$_{{N}}(x={x:.2f})$', fontsize=15)

ax1.set_title(rf'$J_{{n}}(x) = \frac{{2(n-1)}}{{x}} \; J_{{n-1}}(x) - J_{{n-2}}(x)$; x={x:.2f}', fontsize=17, pad=20)
ax2.set_title(rf'$J_{{n-1}}(x) = \frac{{2n}}{{x}} \; J_{{n}}(x) - J_{{n+1}}(x)$; x={x:.2f}', fontsize=17, pad=20)

ax1.set_xticks(n_array); ax2.set_xticks(n_array)
ax3.set_xticks(n_array); ax4.set_xticks(n_array)
ax1.legend(fontsize=15, edgecolor='k')
ax2.legend(fontsize=15, edgecolor='k')
ax3.legend(fontsize=15, edgecolor='k')
ax4.legend(fontsize=15, edgecolor='k')

fig.subplots_adjust(
    top=0.930,
    bottom=0.065,
    left=0.055,
    right=0.985,
    hspace=0.2,
    wspace=0.145
)
plt.show()