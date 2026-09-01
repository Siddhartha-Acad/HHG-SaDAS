! parameters.f90
!
! Fortran counterpart of parameters.py. Mirrors its sections and naming as
! closely as Fortran allows, so any downstream program (GPSM-DSYEV.f90,
! S-matrix-DSYEV.f90, vector_time_evolution.f90, ...) that does `use parameters`
! fetches the SAME system definition (atom, confinement, laser, grids, ...)
! from a single place, just like every *.py script in the project imports
! from parameters.py.
!
! [NOTE] All identifiers that already existed here (n_qn, l_qn, m_qn, N, L,
! l_max, kmax, total_states, Lmap, r_max, r0, alpha_map, dt, check_n_states,
! and the ANSI color strings) are kept byte-for-byte so nothing that already
! `use`s this module needs to change.
!
! [NOTE] The physical conversion constants in the "LASER" section below
! (Int_0_au, hartree_wavelength_nm) are standard atomic-unit constants. If
! your project has its own Atomic_units.f90 (the Fortran analogue of
! Atomic_units.py, which parameters.py imports Int_0/omega_au/T0 from),
! swap these two lines for calls into that module so the numbers match
! bit-for-bit with the Python side.

module parameters
    implicit none

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !           Atom, SAE and Confinement            |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer, parameter :: n_qn = 1, l_qn = 0, m_qn = 0         ! defines initial state. (qn = quantum number)

    character(len=*), parameter :: evolving_atom = 'H'         ! Atoms listed in the 'SAE dataset' section below.
    character(len=*), parameter :: SAE_model = 'SAE-M1'        ! Single active electron model: 'SAE-M1' or 'SAE-M2'.
                                                               ! [NOTE]: For 'Xe' always use 'SAE-M1' (matches parameters.py).

    logical, parameter :: confined = .false.                   ! whether the atom is confined or not?
    character(len=*), parameter :: conf_model = 'P-Gau'        ! confinement potential type, if confined = .true.
                                                               ! options: 'ASW', 'GASW', 'Lor', 'SSW', 'Gau', 'P-Gau'
    logical, parameter :: save_Egvals_with_Smatrix = .true.    ! save eigenvalues alongside the S-matrix (energy level diagram)

    ! 'Confined_atom' / 'Free_atom' sub-directory selector, same role as
    ! Python's `data_dir = 'Confined_atom' if confined else 'Free_atom'`.
    character(len=13), parameter :: data_dir = merge('Confined_atom', 'Free_atom    ', confined)


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                GPSM Parameters                 |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer, parameter :: N = 200            ! P'_N(xj) = 0 ; radial grid size: len(colloc_pt) = N-1
    integer, parameter :: L = 20             ! must be >= l ; angular grid size: len(theta_k) = L+1
    integer, parameter :: l_max = 20         ! Number of partial waves = number of S-matrices = l_max+1
    integer, parameter :: kmax = 50          ! number of GPSM states (maximum k index) in S matrix
    integer, parameter :: total_states = 5   ! how many states you want to keep in the GPSM_state file (.dat)
    real(kind=8), parameter :: Lmap = 80.0d0, r_max = 200.0d0, r0 = 150.0d0    ! radial mapping parameters and absorber radius
    real(kind=8), parameter :: alpha_map = 2.0d0 * Lmap / r_max


    real(kind=8), parameter :: dt = 0.1
    integer, parameter :: check_n_states = 6    ! how many states you want to include in check_Split_operator.f90?


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !          LASER and temporal grid info          |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    real(kind=8), parameter :: lambda_nm = 1064.0d0                  ! wavelength (nm)
    real(kind=8), parameter :: I0 = 5.0d13                           ! Intensity (W/cm^2)

    ! Standard atomic-unit conversion constants (see [NOTE] at top of file):
    real(kind=8), parameter :: Int_0_au = 3.50944758d16              ! atomic unit of intensity (W/cm^2)
    real(kind=8), parameter :: hartree_wavelength_nm = 45.56335d0    ! wavelength (nm) of a 1-Hartree photon

    real(kind=8), parameter :: I0_au = I0 / Int_0_au                 ! Intensity (a.u.)
    real(kind=8), parameter :: E0_au = sqrt(I0_au)                   ! Field amplitude (a.u.)
    real(kind=8), parameter :: w0 = hartree_wavelength_nm / lambda_nm ! Angular frequency (a.u.)
    real(kind=8), parameter :: pi_au = acos(-1.0d0)
    real(kind=8), parameter :: T0 = 2.0d0 * pi_au / w0               ! Period (a.u.)

    integer, parameter      :: cpp = 60                              ! cycles per pulse
    real(kind=8), parameter :: tf = dble(cpp) * T0                   ! total pulse duration (a.u.)

    ! Maximum number of time steps available on the grid t = 0, dt, 2dt, ..., tf
    ! (Fortran analogue of Python's `time_step = len(t) - 1`).
    integer, parameter :: time_step_max = nint(tf / dt)


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !   Time evolution controls: vector_time_evolution.f90   |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    real(kind=8), parameter :: eta_t = 0.00154d0     ! Execution time for a single time-step (dt) evolution (seconds, per-run estimate)
    integer, parameter :: time_step = time_step_max  ! number of time steps desired for evolution. Max possible = time_step_max
    logical, parameter :: show_E_field = .false.     ! whether to display the laser E(t) before evolution starts
    logical, parameter :: print_serial_prog = .true. ! print progress while vector_time_evolution.f90 runs
    integer, parameter :: p_step = 10                ! print every p_step (%) completion


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                  ANSI colors                   |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    character(len=*), parameter :: green  = char(27)//'[1;32m', &
                                   red    = char(27)//'[1;31m', &
                                   yellow = char(27)//'[1;33m', &
                                   white  = char(27)//'[1;37m', &
                                   reset  = char(27)//'[0m'


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                  SAE dataset                   |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Fortran analogue of Python's atomic_params_SAE_M1 / atomic_params_SAE_M2
    ! dicts. Look an atom up with get_SAE_M1('Ar') / get_SAE_M2('Ar').

    type :: sae_m1_t
        character(len=2) :: atom
        real(kind=8) :: Zc, a1, a2, a3, a4, a5, a6
    end type sae_m1_t

    ! Ref: X. M. Tong and C. D. Lin, J. Phys. B: At. Mol. Opt. Phys., 38, 2593 (2005).
    type(sae_m1_t), parameter :: atomic_params_SAE_M1(5) = (/ &
        sae_m1_t('H ', 1.0d0,  0.000d0,  0.000d0,   0.000d0,  0.000d0,  0.000d0, 0.000d0), &
        sae_m1_t('He', 1.0d0,  1.231d0,  0.662d0,  -1.325d0,  1.236d0, -0.231d0, 0.480d0), &
        sae_m1_t('Ne', 1.0d0,  8.069d0,  2.148d0,  -3.570d0,  1.986d0,  0.931d0, 0.602d0), &
        sae_m1_t('Ar', 1.0d0, 16.039d0,  2.007d0, -25.543d0,  4.525d0,  0.961d0, 0.443d0), &
        sae_m1_t('Xe', 1.0d0, 51.356d0,  2.112d0, -99.927d0,  3.737d0,  1.644d0, 0.431d0) /)

    type :: sae_m2_t
        character(len=2) :: atom
        real(kind=8) :: C0, Zc, c, a1, a2, a3, b1, b2, b3
    end type sae_m2_t

    ! Ref: R. Reiff, T. Joyce, A. Jaroń-Becker, and A. Becker, J. Phys. Commun., 4, 065011 (2020).
    type(sae_m2_t), parameter :: atomic_params_SAE_M2(4) = (/ &
        sae_m2_t('H ', 1.0d0,  0.0d0, 0.0000d0,  0.0000d0,   0.0000d0,  0.0000d0, 0.0000d0, 0.0000d0,  0.0000d0), &
        sae_m2_t('He', 1.0d0,  1.0d0, 2.0329d0,  0.3953d0,   0.0000d0,  0.0000d0, 6.1805d0, 0.0000d0,  0.0000d0), &
        sae_m2_t('Ne', 1.0d0,  9.0d0, 0.8870d0, -9.9286d0,  -5.9950d0,  0.0000d0, 1.3746d0, 3.7963d0,  0.0000d0), &
        sae_m2_t('Ar', 1.0d0, 17.0d0, 0.8103d0, -15.9583d0, -27.7467d0, 2.1768d0, 1.2305d0, 4.3946d0, 86.7179d0) /)

contains

    ! Fortran analogue of parameters.py::state_name(n_val, l_val) -> e.g. (1,0) -> "1s"
    function state_symbol(n_val, l_val) result(s)
        integer, intent(in) :: n_val, l_val
        character(len=4) :: s
        character(len=1), parameter :: letters(0:6) = (/ 's', 'p', 'd', 'f', 'g', 'h', 'i' /)
        character(len=1) :: orb

        if (l_val >= 0 .and. l_val <= 6) then
            orb = letters(l_val)
        else
            orb = '?'
        end if
        write(s, '(I0,A1)') n_val, orb
    end function state_symbol

    ! Look up SAE-M1 coefficients by atom symbol, e.g. get_SAE_M1('Ar').
    ! Fortran analogue of Python's atomic_params_SAE_M1[evolving_atom].
    function get_SAE_M1(atom_sym) result(p)
        character(len=*), intent(in) :: atom_sym
        type(sae_m1_t) :: p
        integer :: i

        do i = 1, size(atomic_params_SAE_M1)
            if (trim(atomic_params_SAE_M1(i)%atom) == trim(atom_sym)) then
                p = atomic_params_SAE_M1(i)
                return
            end if
        end do
        print *, 'ERROR: atom "', trim(atom_sym), '" not found in atomic_params_SAE_M1'
        stop 1
    end function get_SAE_M1

    ! Look up SAE-M2 coefficients by atom symbol, e.g. get_SAE_M2('Ar').
    function get_SAE_M2(atom_sym) result(p)
        character(len=*), intent(in) :: atom_sym
        type(sae_m2_t) :: p
        integer :: i

        do i = 1, size(atomic_params_SAE_M2)
            if (trim(atomic_params_SAE_M2(i)%atom) == trim(atom_sym)) then
                p = atomic_params_SAE_M2(i)
                return
            end if
        end do
        print *, 'ERROR: atom "', trim(atom_sym), '" not found in atomic_params_SAE_M2'
        stop 1
    end function get_SAE_M2

end module parameters