program read_GPSM_data
    implicit none

    integer :: i
    integer, parameter :: N = 200, kmax = 5
    real(kind=8), dimension(N-1) :: rx
    real(kind=8), dimension(N-1, kmax) :: A

    open(unit=12, file='../GPSM-DSYEV_states.bin', form='unformatted', access='stream', status='old')

    do i = 1, N-1
        read(12) rx(i), A(i, 1:kmax)
    end do

    close(12)

    print *, sum(rx * A(:, 1) ** 2)

end program read_GPSM_data

