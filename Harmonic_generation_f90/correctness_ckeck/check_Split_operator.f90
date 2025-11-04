program check_Split_operator
    use parameters
    implicit none

    complex*16 :: S_matrix(N-1, N-1)
    complex*16 :: S_all(0:l_max, N-1, N-1)
    integer :: l_ind

    open(unit=21, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices.bin', &
        form='unformatted', access='sequential', status='old')

    do l_ind = 0, l_max
        read(21) S_matrix
        S_all(l_ind, :, :) = S_matrix
        print '(A, I3, A)', 'Loaded S(l=', l_ind + m_qn, ')'
    end do

    close(21)

end program check_Split_operator
