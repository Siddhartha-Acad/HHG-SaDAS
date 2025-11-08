program check_Split_operator
    use parameters
    implicit none

    integer :: i, recl_size
    integer, parameter :: l_ind = 1             ! select S_matrix, compatible with GPSM l_qn.
    complex(kind=8) :: S_matrix(N-1, N-1)

    real(kind=8) :: max_abs_err
    real(kind=8), dimension(N-1) :: r, A_tilde
    real(kind=8), dimension(N-1, total_states) :: A

    inquire(iolength=recl_size) S_matrix            ! "How many units does S_matrix need?"

    open(unit=10, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
            form='unformatted', access='stream', status='old')
        read(10) r, A
    close(10)

    open(unit=11, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
            form='unformatted', access='direct', recl=recl_size, status='old')
        read(11, rec=l_ind) S_matrix        ! rec = 1, 2, 3, ..., l_max
    close(11)

    print '(A11, A, A11, A, A6)', 'S(l) @ A(n)', ' | ', 'abs_max_err', ' | ', 'status'
    print '(A)', repeat('-', 35)

    do i = 1, total_states
        A_tilde = matmul(S_matrix, cmplx(A(:, i), 0.0d0, kind=8))
        max_abs_err = maxval(abs(A_tilde - A(:, i)))

        if (max_abs_err .lt. 1.0d-3) then
            print '(A,I0,A,I0,A,A,ES11.3,A,A)', 'S(', l_ind-1, ') @ A(', i, ')', ' | ', max_abs_err, ' | ', '[PASS]'
        else
            print '(A,I0,A,I0,A,A,ES11.3,A,A)', 'S(', l_ind-1, ') @ A(', i, ')', ' | ', max_abs_err, ' | ', '[FAIL]'
        end if
    end do

end program check_Split_operator
