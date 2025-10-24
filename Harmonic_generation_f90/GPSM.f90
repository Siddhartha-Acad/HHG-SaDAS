! fortran implementation for solving 
! time independent Schrodinger equation using GPSM

program GPSM
    implicit none
    integer :: i
    integer, parameter :: N = 200
    real(kind=8), dimension(N-1) :: x   ! collocation points
    real(kind=8), external :: f_p

    !~~~~~~~~~: reading collocation points :~~~~~~~~~
    open(unit=10, file='./collocation_points/generator/Algo-3_N=200_Gauss_Lobatto_collocation_points.dat', &
         status='old', action='read')

    do i=1, N-1
        read(10, *) x(i)
    end do
    close(10)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, f_p(x(1))
    

end program GPSM


pure real(kind=8) function f_p(x)
	implicit none
	real(kind=8), intent(in) :: x
	real(kind=8) :: r_max, Lmap, alpha_map
	
	r_max = 200; Lmap = 80
	alpha_map = 2 * Lmap / r_max
    f_p = Lmap * (alpha_map + 2) / (1 - x + alpha_map)**2
end function f_p
