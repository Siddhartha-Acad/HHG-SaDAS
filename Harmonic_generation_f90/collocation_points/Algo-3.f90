! Newton-Raphson method of root finding.

program main
    implicit none
    real (kind=8) :: x_i, root

    x_i = 2.5d0

    call newton_raphson(x_i, root, .false.)
    print *, 'root = ', root
end program main


subroutine newton_raphson(x_i, root, debug)
    implicit none
    real (kind=8), intent(in) :: x_i        ! initial guess value
    real (kind=8), intent(out) :: root
    logical, intent(in) :: debug
    real (kind=8) :: x_old, x_new
    real (kind=8) :: tol, rtol
    integer :: n

    tol = 1.0d-14       ! absolute error tolerance
    rtol = 0.0d0        ! relative tolerance
    x_old = x_i

    if (debug) then
        print '(A3, 2A20, A25)', &
            'n', 'x_n', 'x_{n+1}', 'err = |x_{n+1} - x_n|'
    end if

    do n = 1, 100
        x_new = x_old - sin(x_old) / cos(x_old)

        if (debug) then
            print '(I3, 2F20.16, E25.16)', &
                n, x_old, x_new, abs(x_new - x_old)
        end if

        if (abs(x_new - x_old) .lt. (tol + rtol*abs(x_new))) then
            root = x_new
            return
        end if

        x_old = x_new
    end do

    ! Reaching here is signature that convergence failed.
    root = 1.0d2
end subroutine newton_raphson
