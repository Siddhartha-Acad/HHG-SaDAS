program read_GPSM_data
    implicit none
    integer :: i
    integer, parameter :: N = 200, kmax = 5, l = 0
    real :: rad_theory
    real(kind=8), dimension(N-1) :: r
    real(kind=8), dimension(N-1, kmax) :: A

    open(unit=12, file='../GPSM-DSYEV_states.bin', form='unformatted', access='stream', status='old')
        do i = 1, N-1
            read(12) r(i), A(i, 1:kmax)
        end do
    close(12)
    
    print '(A2, A14, A22)', 'n', 'radius_Th', 'radius_GPSM'
    do i = 1, kmax
        rad_theory = 0.5 * real(3*i**2 - l*(l+1))
        print '(I2, F14.8, F22.16)', i, rad_theory, sum(r * A(:, i)**2)
    end do

end program read_GPSM_data

