"""
HHG spectra analysis for Fortran time-evolution data.

Reads:
    t, E(t), d(t), Ps(t)

and plots:
    1. Dipole moment d(t)
    2. HHG spectrum
    3. Electric field E(t) and survival probability Ps(t)
    4. Survival probability Ps(t)
    5. Ionisation probability Pi(t)
"""

import numpy as np
from pathlib import Path
import scipy.fft as ft
import matplotlib.pyplot as plt


# ============================================================
#                 Laser parameters
# ============================================================

# Set this to the fundamental laser frequency used
# in your Fortran calculation.
w0 = 0.057       # a.u.  <-- CHANGE if required

T = 2 * np.pi / w0


# ============================================================
#                 Figure settings
# ============================================================

width = 6.2
height = 4
fig_scale_factor = 2
fig_size = (fig_scale_factor * width,
            fig_scale_factor * height)

plt.rc('font', family='serif', size=14)


# ============================================================
#                 Load evolution data
# ============================================================

def load_evo_data(path):

    path = Path(path)

    data = np.loadtxt(path, skiprows=1)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    return data


# ============================================================
#                         Main
# ============================================================

if __name__ == '__main__':

    evo_dir = Path(__file__).resolve().parent.parent / 'Time_evolution_data'

    files = sorted(evo_dir.glob('*.dat'))

    if not files:
        print('No time-evolution .dat files found in ./Time_evolution_data')
        raise SystemExit(0)

    # Use the first file
    evo_file = files[0]

    data = load_evo_data(evo_file)

    print(f'Loaded: {evo_file.name}')
    print(f'shape = {data.shape}')

    # --------------------------------------------------------
    # Extract columns
    # --------------------------------------------------------

    t = data[:, 0]
    E_t = data[:, 1]
    dipole_moment = data[:, 2]
    survival_probability = data[:, 3]

    # Determine dt directly from the data
    dt = np.mean(np.diff(t))

    # Total evolution time
    tf = t[-1] - t[0]

    print(f'Time step dt       = {dt:.6f} a.u.')
    print(f'Total time         = {tf:.6f} a.u.')
    print(f'Number of points   = {len(t)}')
    print(f'Fundamental w0     = {w0:.6f} a.u.')
    print(f'Optical period T   = {T:.6f} a.u.')


    # ========================================================
    #                    HHG spectrum
    # ========================================================

    dipole_mom_FFT = ft.fft(dipole_moment)

    freq_w = ft.fftfreq(
        len(t),
        d=dt
    ) * 2 * np.pi

    dip_mom_power_spectra = np.abs(
        dipole_mom_FFT / tf
    )**2

    # Positive frequencies only
    positive_freq_mask = freq_w > 0

    freq_pos = freq_w[positive_freq_mask]

    dip_mom_power_spectra = \
        dip_mom_power_spectra[positive_freq_mask]

    # Harmonic order
    harmonic_order = freq_pos / w0


    # ========================================================
    #                Conversion efficiency
    # ========================================================

    harmonic_energy = np.trapezoid(
        dip_mom_power_spectra,
        freq_pos
    )

    laser_energy = np.trapezoid(
        E_t**2,
        t
    )

    eta = harmonic_energy / laser_energy


    print("\n~~~~~~~~~~~: Harmonic Spectra Summary :~~~~~~~~~~")
    print(f"Number of time points : {len(t)}")
    print(f"Time step dt          : {dt:.6f} a.u.")
    print(f"Fundamental frequency : {w0:.6f} a.u.")
    print(f"Nyquist frequency     : {freq_pos[-1]:.6f} a.u.")

    print("\n~~~~~~~~~~~~: Conversion efficiency :~~~~~~~~~~~~")
    print(f"Total laser energy    : {laser_energy:.6e}")
    print(f"Total harmonic energy : {harmonic_energy:.6e}")
    print(f"Conversion efficiency : {eta:.6e}")


    # ========================================================
    #                         Figures
    # ========================================================

    fig1 = plt.figure(figsize=fig_size)
    fig2 = plt.figure(figsize=fig_size)

    # --------------------------------------------------------
    # Figure 1
    # --------------------------------------------------------

    ax1 = fig1.add_subplot(211)
    ax2 = fig1.add_subplot(212)

    # --------------------------------------------------------
    # Figure 2
    # --------------------------------------------------------

    ax3 = fig2.add_subplot(211)
    ax4 = ax3.twinx()

    ax5 = fig2.add_subplot(223)
    ax6 = fig2.add_subplot(224)


    # ========================================================
    #                         Plots
    # ========================================================

    # Dipole moment
    ax1.plot(
        t / T,
        dipole_moment,
        lw=2,
        label='d(t)'
    )

    # HHG spectrum
    ax2.plot(
        harmonic_order,
        np.log10(dip_mom_power_spectra),
        label=r'log$_{10}$[P($\omega$)]'
    )

    # Electric field
    ax3.plot(
        t / T,
        E_t,
        label='E(t)'
    )

    # Survival probability on secondary axis
    ax4.plot(
        t / T,
        survival_probability,
        label=r'P$_s$(t)'
    )

    # Survival probability
    ax5.plot(
        t / T,
        survival_probability,
        lw=2,
        label=r'P$_s$(t)'
    )

    # Ionisation probability
    ionisation_probability = 1 - survival_probability

    ax6.plot(
        t / T,
        ionisation_probability,
        lw=2,
        label=r'P$_i$(t)'
    )


    # ========================================================
    #                         Labels
    # ========================================================

    ax1.set_xlabel('Optical cycles (t/T)')
    ax1.set_ylabel('d(t) (a.u.)')

    ax2.set_xlabel(r'Harmonic Order ($\omega/\omega_0$)')
    ax2.set_ylabel(r'log$_{10}$[P($\omega$)]')

    ax3.set_xlabel('Optical cycles (t/T)')
    ax3.set_ylabel('E(t) (a.u.)')

    ax4.set_ylabel(r'P$_s$(t)')

    ax5.set_xlabel('Optical cycles (t/T)')
    ax5.set_ylabel(r'P$_s$(t)')

    ax6.set_xlabel('Optical cycles (t/T)')
    ax6.set_ylabel(r'P$_i$(t)')


    # ========================================================
    #                         Legends
    # ========================================================

    ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    ax3.legend(loc='upper right')
    ax4.legend(loc='lower left')
    ax5.legend(loc='upper right')
    ax6.legend(loc='lower right')


    # ========================================================
    #                         Limits
    # ========================================================

    ax2.set_xlim(-1, min(102, harmonic_order[-1]))

    # Avoid log10(0)
    ax2.set_ylim(
        np.nanmin(np.log10(dip_mom_power_spectra)),
        np.nanmax(np.log10(dip_mom_power_spectra))
    )

    ax5.set_ylim(0, 1.05)
    ax6.set_ylim(0, 1.05)


    # ========================================================
    #                         Titles
    # ========================================================

    fig1.suptitle(evo_file.name, fontsize=12)
    fig2.suptitle(evo_file.name, fontsize=12)

    fig1.subplots_adjust(
        top=0.94,
        bottom=0.069,
        right=0.98,
        left=0.079,
        hspace=0.195
    )

    fig2.subplots_adjust(
        top=0.94,
        bottom=0.069,
        right=0.92,
        left=0.08,
        hspace=0.195,
        wspace=0.074
    )


    plt.show()