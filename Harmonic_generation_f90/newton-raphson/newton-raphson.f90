! Newton-Raphson method of root finding.

program main
    implicit none
    integer :: n
    real (kind=8) :: x_i, x_n      ! init_guess, iteration_val

    x_i = 3.0d0

    do n = 1, 10
        x_n = x_i - sin(x_i) / cos(x_i)

        x_i = x_n
    end do

    print *, x_i, 



end program main