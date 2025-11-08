program check_Split_operator
    use parameters
    implicit none

    integer :: i, recl_size
    integer, parameter :: l_ind = 1             ! select S_matrix, compatible with GPSM l_qn.
    complex(kind=8) :: S_matrix(N-1, N-1)

    real(kind=8) :: rel_err
    real(kind=8), dimension(N-1) :: r, A_tilde
    real(kind=8), dimension(N-1, total_states) :: A

    inquire(iolength=recl_size) S_matrix            ! "How many units does S_matrix need?"
    print *, "S-Matrix Record length:", recl_size

    open(unit=10, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
            form='unformatted', access='stream', status='old')
        read(10) r, A
    close(10)

    open(unit=11, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
            form='unformatted', access='direct', recl=recl_size, status='old')
        read(11, rec=l_ind) S_matrix        ! rec = 1, 2, 3, ..., l_max
    close(11)

    A_tilde = matmul(S_matrix, cmplx(A(:, 1), 0.0d0, kind=8))
    rel_err = sum(abs(A_tilde**2 - A(:, 1)**2) / abs(A(:, 1)**2)) / (N-1)

    print *, ''
    print *, '=== Error Analysis ==='
    print *, 'Max absolute error:     ', maxval(abs(A_tilde - A(:, 1)))
    print *, 'Relative error (norms): ', abs(norm2(A_tilde) - norm2(A(:, 1))) / norm2(A(:, 1))


end program check_Split_operator
