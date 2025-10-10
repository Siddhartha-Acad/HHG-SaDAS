! Newton-Raphson method of root finding.

program main
    implicit none
    real (kind=8) :: x_i
    real (kind=8) :: root, newton_raphson

    x_i = 0.5d0

    root = newton_raphson(x_i)

    print *, 'root = ', root

end program main



function newton_raphson(x_i) result(root)
    implicit none
    real (kind=8), intent(in) :: x_i        ! initial guess value
    real (kind=8) :: root, x_old, x_new
    real (kind=8) :: tol, rtol, rel_err
    integer :: n

    tol = 1.0d-14
    rtol = 0.0d0

    x_old = x_i
    print '(a3, 2a20, a25)', 'n', 'x_old', 'x_new', 'err'

    do n = 1, 100
        x_new = x_old - sin(x_old) / cos(x_old)

        print '(I3, 2F20.16, E25.16)', n, x_old, x_new, abs(x_new - x_old)

        if (abs(x_new - x_old) .lt. (tol + rtol*abs(x_new))) then
            root = x_new
            return
        end if

        x_old = x_new
    end do

    ! Reaching here is signature that convergence failed.
    root = 1.0d2
end function newton_raphson