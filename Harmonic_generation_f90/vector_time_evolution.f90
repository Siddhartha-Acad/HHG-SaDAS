! vector_time_evolution.f90 : Fortran implementation of
! ../Harmonic_generation_py/vector_time_evolution.py
!
! This version follows the split-operator / partial-wave evolution used in the
! Python script: each step applies the S(l) propagation in the radial basis,
! reconstructs the angular wavefunction, multiplies by exp(-i V_int dt), and
! projects back to partial waves before the absorber is applied.

program main
    use parameters
    use timer_mod
    use legendre_stuff, only: a_legendre
    implicit none

    integer :: i, j, l_idx, ti, rec, info
    integer :: S_recl_size, state_recl_size, time_step
    integer(kind=8) :: state_offset
    integer(kind=8), external :: factorial

    real(kind=8), parameter :: pi = acos(-1.0d0)
    real(kind=8), parameter :: E0_au = 0.1d0
    real(kind=8), parameter :: w0 = 0.057d0
    real(kind=8), parameter :: cpp = 60.0d0
    real(kind=8), parameter :: laser_dt = 0.1d0

    real(kind=8), dimension(L+1) :: roots, weights, cos_theta
    real(kind=8), dimension(N-1) :: x_glob, r, absorber, A_r
    real(kind=8), dimension(N-1) :: dipole_vals, population_vals
    real(kind=8), dimension(L+1, N-1) :: V_int_matrix
    real(kind=8), dimension(L+1, N-1) :: psi_1, psi_2
    real(kind=8), dimension(N-1, total_states) :: state_block
    real(kind=8), dimension(N-1, N-1) :: S_l_matrix
    complex(kind=8), dimension(l_max+1, N-1) :: init_glm, glm_tilde, tmp_glm
    complex(kind=8), dimension(l_max+1, N-1) :: A_mat
    complex(kind=8), dimension(L+1, N-1) :: psi_tmp
    complex(kind=8), dimension(l_max+1, N-1, N-1) :: S_matrix_all
    complex(kind=8), dimension(N-1) :: matvec

    real(kind=8) :: exec_time, t_mid, E_val
    integer :: n_steps


    ! ---------------------------------------------------------------
    ! 1) Read Gauss-Legendre angular nodes/weights for the theta grid
    ! ---------------------------------------------------------------
    open(unit=20, file='../collocation_points/generator/Gauss_Legendre_collocation_points_and_weights.dat', &
         status='old', action='read')
        read(20, *) roots, weights
    close(20)

    cos_theta = roots

    ! ---------------------------------------------------------------
    ! 2) Read Gauss-Lobatto radial grid points from the collocation data
    ! ---------------------------------------------------------------
    open(unit=21, file='../collocation_points/generator/Algo-3_Gauss_Lobatto_collocation_points.dat', &
         status='old', action='read')
        read(21, *) x_glob
    close(21)

    r = Lmap * (1.0d0 + x_glob) / (1.0d0 - x_glob + alpha_map)

    ! ---------------------------------------------------------------
    ! 3) Read the initial GPSM eigenvector from the generated stream file
    !    The file layout is: f(r), H(:,1:total_states)
    ! ---------------------------------------------------------------
    open(unit=22, file='GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
         form='unformatted', access='stream', status='old')
        read(22) r
        read(22) state_block
    close(22)

    A_r = state_block(:, 1)
    init_glm = (0.0d0, 0.0d0)
    init_glm(1, :) = cmplx(A_r, 0.0d0, kind=8)

    ! ---------------------------------------------------------------
    ! 4) Read all S(l) matrices from the direct-access binary file
    ! ---------------------------------------------------------------
    inquire(iolength=S_recl_size) S_l_matrix
    open(unit=23, file='GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
         form='unformatted', access='direct', recl=S_recl_size, status='old')

    do rec = 1, l_max + 1
        read(23, rec=rec) S_l_matrix
        S_matrix_all(rec, :, :) = cmplx(S_l_matrix, 0.0d0, kind=8)
    end do
    close(23)

    ! ---------------------------------------------------------------
    ! 5) Build the absorber mask and a small time grid for propagation
    ! ---------------------------------------------------------------
    do i = 1, N-1
        if (0.0d0 < r(i) .and. r(i) <= r0) then
            absorber(i) = 1.0d0
        else if (r(i) > r0 .and. r(i) < r_max) then
            absorber(i) = cos(pi * (r(i) - r0) / (2.0d0 * (r_max - r0)))**0.25d0
        else
            absorber(i) = 0.0d0
        end if
    end do

    ! For demonstration and stable testing, use a short propagation window.
    ! For the full production run you can replace this with the Python value:
    !   time_step = int((cpp * 2*pi / w0) / laser_dt)
    time_step = 10

    call tick()

    ! ---------------------------------------------------------------
    ! 6) Main split-operator time evolution loop
    ! ---------------------------------------------------------------
    do ti = 1, time_step
        ! 6a) Apply the radial propagation operator in each l channel
        do l_idx = 1, l_max + 1
            matvec = matmul(S_matrix_all(l_idx, :, :), init_glm(l_idx, :))
            A_mat(l_idx, :) = matvec
        end do

        ! 6b) Reconstruct psi(theta, r) = sum_l Y_lm(theta) * A_l(r)
        psi_tmp = (0.0d0, 0.0d0)
        do j = 1, L+1
            do i = 1, N-1
                do l_idx = 1, l_max + 1
                    psi_tmp(j, i) = psi_tmp(j, i) + cmplx(a_legendre(l_idx-1, 0, roots(j)), 0.0d0, kind=8) * A_mat(l_idx, i)
                end do
            end do
        end do

        ! 6c) Interaction step: exp(-i V_int dt)
        t_mid = (dble(ti) - 0.5d0) * laser_dt
        E_val = E0_au * sin(w0 * t_mid) * sin(w0 * t_mid / (2.0d0 * cpp))**2

        do j = 1, L+1
            do i = 1, N-1
                V_int_matrix(j, i) = -E_val * cos_theta(j) * r(i)
                psi_2(j, i) = exp(cmplx(0.0d0, -V_int_matrix(j, i) * laser_dt, kind=8)) * psi_tmp(j, i)
            end do
        end do

        ! 6d) Project back to partial waves g_lm(r)
        glm_tilde = (0.0d0, 0.0d0)
        do l_idx = 1, l_max + 1
            do i = 1, N-1
                tmp_glm(l_idx, i) = (0.0d0, 0.0d0)
                do j = 1, L+1
                    tmp_glm(l_idx, i) = tmp_glm(l_idx, i) + &
                        cmplx(weights(j) * a_legendre(l_idx-1, 0, roots(j)), 0.0d0, kind=8) * psi_2(j, i)
                end do
                glm_tilde(l_idx, i) = tmp_glm(l_idx, i) / cmplx(norm_coeff(l_idx-1, 0), 0.0d0, kind=8)
            end do
        end do

        ! 6e) Propagate one more half-step and apply absorber
        do l_idx = 1, l_max + 1
            matvec = matmul(S_matrix_all(l_idx, :, :), glm_tilde(l_idx, :))
            init_glm(l_idx, :) = matvec * absorber
        end do

        ! 6f) Diagnostics for this step
        dipole_vals(ti) = dipole_scalar(r, init_glm)
        population_vals(ti) = population_scalar(init_glm)
    end do

    call tock(exec_time)

    print *, 'Vector time-evolution finished.'
    print '(A,F12.6,A)', 'Wall time : ', exec_time, ' s'
    print '(A,I0)', 'Time steps used : ', time_step
    print '(A,F12.6)', 'Final dipole moment : ', dipole_vals(time_step)
    print '(A,F12.6)', 'Final population    : ', population_vals(time_step)

contains

    real(kind=8) function norm_coeff(l_val, m_val)
        integer, intent(in) :: l_val, m_val
        integer(kind=8) :: fact_n, fact_p
        integer(kind=8), external :: factorial
        real(kind=8) :: phase, pi_local
        pi_local = acos(-1.0d0)

        if (l_val < 0 .or. m_val < 0 .or. m_val > l_val) then
            norm_coeff = 0.0d0
            return
        end if

        fact_n = factorial(l_val - m_val)
        fact_p = factorial(l_val + m_val)
        phase = merge(1.0d0, -1.0d0, mod(abs(m_val), 2) == 0)
        norm_coeff = phase * sqrt((2.0d0 * l_val + 1.0d0) * dble(fact_n) / (4.0d0 * pi_local * dble(fact_p)))
    end function norm_coeff

    real(kind=8) function dipole_scalar(rp, glm_arr)
        real(kind=8), intent(in) :: rp(N-1)
        complex(kind=8), intent(in) :: glm_arr(l_max+1, N-1)
        integer :: l
        real(kind=8) :: sum_r, alpha_l

        sum_r = 0.0d0
        do l = 1, l_max
            alpha_l = sqrt(((dble(l))**2) / ((2.0d0 * dble(l) - 1.0d0) * (2.0d0 * dble(l) + 1.0d0)))
            sum_r = sum_r + alpha_l * sum(rp(1:N-1) * real(conjg(glm_arr(l, :)) * glm_arr(l+1, :)))
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

    real(kind=8) function E_field_value(tval)
        real(kind=8), intent(in) :: tval
        E_field_value = E0_au * sin(w0 * tval) * sin(w0 * tval / (2.0d0 * cpp))**2
    end function E_field_value

end program main
