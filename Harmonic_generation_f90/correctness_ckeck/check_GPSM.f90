program read_GPSM_data
    use parameters
    implicit none
    integer :: i
    real :: rad_theory
    real(kind=8), dimension(N-1) :: r
    real(kind=8), dimension(N-1, kmax) :: A

    open(unit=12, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
        form='unformatted', access='stream', status='old')
        read(12) r, A
    close(12)

    print '(A2, A14, A22)', 'n', 'radius_Th', 'radius_GPSM'
    do i = 1, kmax
        rad_theory = 0.5 * real(3*i**2 - l_qn*(l_qn+1))
        print '(I2, F14.8, F22.16)', i, rad_theory, sum(r * A(:, i)**2)
    end do

end program read_GPSM_data

