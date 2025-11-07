program check_Split_operator
    use parameters
    implicit none

    complex*16 :: S_matrix(N-1, N-1)
    complex*16 :: S_all(0:l_max, N-1, N-1)
    integer :: l_ind

    real(kind=8), dimension(N-1, kmax) :: A

    open(unit=10, file='../GPSM-DSYEV_states.bin', form='unformatted', access='stream', status='old')
        read(10) r, A
    close(10)

    open(unit=11, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices.bin', &
        form='unformatted', access='sequential', status='old')

    do l_ind = 0, l_max
        read(11) S_matrix
        S_all(l_ind, :, :) = S_matrix
        print '(A, I3, A)', 'Loaded S(l=', l_ind + m_qn, ')'
    end do

    close(11)

end program check_Split_operator
