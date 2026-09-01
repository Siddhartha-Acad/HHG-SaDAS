! vector_time_evolution.f90 : Fortran implementation of
! ../Harmonic_generation_py/vector_time_evolution.py
!
! This version reproduces the short validation loop used in the Python
! reference: apply S(l), reconstruct psi(theta,r), multiply by the interaction
! phase, project back to g_lm(r), and absorb the outer region.

program main
    use parameters
    use timer_mod
    use legendre_stuff, only: a_legendre
    implicit none

    integer :: i, j, l_idx, ti, rec
    integer :: S_recl_size, time_step
    integer(kind=8), external :: factorial
    real(kind=8), external :: N_fact, C_fact
    character(len=256) :: colloc_file, state_file, smat_file, gl_file

    real(kind=8), parameter :: pi = acos(-1.0d0)
    real(kind=8), parameter :: E0_au = 0.1d0
    real(kind=8), parameter :: w0 = 0.057d0
    real(kind=8), parameter :: cpp = 60.0d0
    real(kind=8), parameter :: laser_dt = 0.1d0

    real(kind=8), dimension(L+1) :: roots, weights, cos_theta
    real(kind=8), dimension(N-1) :: x_glob, r, absorber, A_r
    real(kind=8), dimension(L+1, l_max) :: Y_T
    real(kind=8), dimension(N-1, total_states) :: state_block
    complex(kind=8), dimension(N-1, N-1) :: S_l_matrix_c
    complex(kind=8), dimension(l_max+1, N-1) :: init_glm, glm_tilde, tmp_glm, A_mat
    complex(kind=8), dimension(L+1, N-1) :: psi_1, psi_2
    complex(kind=8), dimension(l_max+1, N-1, N-1) :: S_matrix_all
    complex(kind=8), dimension(N-1) :: matvec

    real(kind=8) :: exec_time, t_mid, E_val, norm_factor
    real(kind=8), allocatable :: dipole_vals(:), population_vals(:)

    ! Short validation run, matching the Python script's fixed 10-step check.
    time_step = 10
    allocate(dipole_vals(time_step), population_vals(time_step))

    ! ---------------------------------------------------------------
    ! 1) Radial grid: N-1 Gauss-Lobatto collocation points
    ! ---------------------------------------------------------------
    colloc_file = 'collocation_points/generator/Algo-3_Gauss_Lobatto_collocation_points.dat'
    if (.not. file_exists(colloc_file)) then
        colloc_file = '../collocation_points/generator/Algo-3_Gauss_Lobatto_collocation_points.dat'
    end if
    if (.not. file_exists(colloc_file)) then
        print *, 'ERROR: collocation file not found. Checked:'
        print *, '  - ', trim(colloc_file)
        stop 1
    end if

    open(unit=21, file=trim(colloc_file), status='old', action='read')
        read(21, *) x_glob
    close(21)

    ! ---------------------------------------------------------------
    ! 2) Angular grid: Gauss-Legendre nodes/weights for the theta grid
    ! ---------------------------------------------------------------
    gl_file = 'collocation_points/generator/Gauss_Legendre_collocation_points_and_weights.dat'
    if (.not. file_exists(gl_file)) then
        gl_file = '../collocation_points/generator/Gauss_Legendre_collocation_points_and_weights.dat'
    end if
    if (.not. file_exists(gl_file)) then
        print *, 'ERROR: Gauss-Legendre collocation/weight file not found. Checked:'
        print *, '  - ', trim(gl_file)
        stop 1
    end if

    open(unit=20, file=trim(gl_file), status='old', action='read')
        read(20, *) roots, weights
    close(20)

    r = Lmap * (1.0d0 + x_glob) / (1.0d0 - x_glob + alpha_map)
    cos_theta = roots

    state_file = 'GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin'
    if (.not. file_exists(state_file)) then
        state_file = '../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin'
    end if
    if (.not. file_exists(state_file)) then
        print *, 'ERROR: GPSM state file not found. Checked:'
        print *, '  - ', trim(state_file)
        stop 1
    end if

    open(unit=22, file=trim(state_file), form='unformatted', access='stream', status='old')
        read(22) r
        read(22) state_block
    close(22)

    A_r = state_block(:, 1)
    init_glm = (0.0d0, 0.0d0)
    init_glm(1, :) = cmplx(A_r, 0.0d0, kind=8)

    inquire(iolength=S_recl_size) S_l_matrix_c
    smat_file = 'GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin'
    if (.not. file_exists(smat_file)) then
        smat_file = '../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin'
    end if
    if (.not. file_exists(smat_file)) then
        print *, 'ERROR: S-matrix file not found. Checked:'
        print *, '  - ', trim(smat_file)
        stop 1
    end if

    open(unit=23, file=trim(smat_file), form='unformatted', access='direct', recl=S_recl_size, status='old')
    do rec = 1, l_max + 1
        read(23, rec=rec) S_l_matrix_c
        S_matrix_all(rec, :, :) = S_l_matrix_c
    end do
    close(23)

    do i = 1, N-1
        if (0.0d0 < r(i) .and. r(i) <= r0) then
            absorber(i) = 1.0d0
        else if (r(i) > r0 .and. r(i) < r_max) then
            absorber(i) = cos(pi * (r(i) - r0) / (2.0d0 * (r_max - r0)))**0.25d0
        else
            absorber(i) = 0.0d0
        end if
    end do

    do j = 1, L+1
        do l_idx = 1, l_max
            Y_T(j, l_idx) = N_fact(l_idx-1, 0) * a_legendre(l_idx-1, 0, roots(j))
        end do
    end do

    call tick()

    do ti = 1, time_step
        do l_idx = 1, l_max
            matvec = matmul(S_matrix_all(l_idx, :, :), init_glm(l_idx, :))
            A_mat(l_idx, :) = matvec
        end do

        psi_1 = (0.0d0, 0.0d0)
        do j = 1, L+1
            do i = 1, N-1
                do l_idx = 1, l_max
                    psi_1(j, i) = psi_1(j, i) + cmplx(Y_T(j, l_idx), 0.0d0, kind=8) * A_mat(l_idx, i)
                end do
            end do
        end do

        t_mid = (dble(ti) - 0.5d0) * laser_dt
        E_val = E0_au * sin(w0 * t_mid) * sin(w0 * t_mid / (2.0d0 * cpp))**2

        do j = 1, L+1
            do i = 1, N-1
                psi_2(j, i) = exp(cmplx(0.0d0, -E_val * cos_theta(j) * r(i) * laser_dt, kind=8)) * psi_1(j, i)
            end do
        end do

        glm_tilde = (0.0d0, 0.0d0)
        do l_idx = 1, l_max
            do i = 1, N-1
                tmp_glm(l_idx, i) = (0.0d0, 0.0d0)
                do j = 1, L+1
                    tmp_glm(l_idx, i) = tmp_glm(l_idx, i) + &
                        cmplx(weights(j) * a_legendre(l_idx-1, 0, roots(j)), 0.0d0, kind=8) * psi_2(j, i)
                end do
                norm_factor = N_fact(l_idx-1, 0) * C_fact(l_idx-1, 0)
                glm_tilde(l_idx, i) = tmp_glm(l_idx, i) / cmplx(norm_factor, 0.0d0, kind=8)
            end do
        end do

        do l_idx = 1, l_max
            matvec = matmul(S_matrix_all(l_idx, :, :), glm_tilde(l_idx, :))
            init_glm(l_idx, :) = matvec * absorber
        end do

        dipole_vals(ti) = dipole_scalar(r, init_glm)
        population_vals(ti) = population_scalar(init_glm)

        print '(A,I3,A,ES16.8)', 'step ', ti, ' | d(t) = ', dipole_vals(ti)
    end do

    call tock(exec_time)

    print *, 'Vector time-evolution finished.'
    print '(A,F12.6,A)', 'Wall time : ', exec_time, ' s'
    print '(A,I0)', 'Time steps used : ', time_step
    print '(A,F12.6)', 'Final dipole moment : ', dipole_vals(time_step)
    print '(A,F12.6)', 'Final population    : ', population_vals(time_step)

contains

    logical function file_exists(path)
        character(len=*), intent(in) :: path
        inquire(file=trim(path), exist=file_exists)
    end function file_exists

    real(kind=8) function dipole_scalar(rp, glm_arr)
        real(kind=8), intent(in) :: rp(N-1)
        complex(kind=8), intent(in) :: glm_arr(l_max+1, N-1)
        integer :: l
        real(kind=8) :: sum_r, alpha_l
        sum_r = 0.0d0
        do l = 1, l_max
            alpha_l = sqrt(((dble(l))**2) / ((2.0d0 * dble(l) - 1.0d0) * (2.0d0 * dble(l) + 1.0d0)))
            sum_r = sum_r + alpha_l * sum(rp * real(conjg(glm_arr(l, :)) * glm_arr(l+1, :)))
        end do
        dipole_scalar = 2.0d0 * sum_r
    end function dipole_scalar

    real(kind=8) function population_scalar(glm_arr)
        complex(kind=8), intent(in) :: glm_arr(l_max+1, N-1)
        integer :: i, l
        population_scalar = 0.0d0
        do l = 1, l_max + 1
            do i = 1, N-1
                population_scalar = population_scalar + abs(glm_arr(l, i))**2
            end do
        end do
    end function population_scalar

end program main