program check_Split_operator
    use parameters
    implicit none

    integer :: l_ind, recl_size
    complex(kind=8) :: S_matrix(N-1, N-1)

    real(kind=8), dimension(N-1) :: r
    real(kind=8), dimension(N-1, total_states) :: A

    inquire(iolength=recl_size) S_matrix            ! "How many units does S_matrix need?"

    open(unit=10, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
        form='unformatted', access='stream', status='old')
        read(10) r, A
    close(10)

    open(unit=11, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
            form='unformatted', access='direct', recl=recl_size, status='old')
        read(11, rec=l_ind) S_matrix
        print '(A, I3, A)', 'Loaded S(l=', l_ind + m_qn - 1, ')'
    close(11)

    
end program check_Split_operator
