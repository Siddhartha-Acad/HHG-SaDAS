! fortran implementation for solving 
! time independent Schrodinger equation using GPSM

program GPSM
    implicit none
    integer :: i
    integer, parameter :: N = 200
    real(kind=8), dimension(N-1) :: x   ! collocation points

    !~~~~~~~~~: reading collocation points :~~~~~~~~~
    open(unit=10, file='./collocation_points/generator/Algo-3_N=200_Gauss_Lobatto_collocation_points.dat', &
         status='old', action='read')

    do i=1, N-1
        read(10, *) x(i)
    end do

    close(10)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, x(199)


end program GPSM
