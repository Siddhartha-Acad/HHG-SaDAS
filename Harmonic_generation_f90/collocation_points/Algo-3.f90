! Newton-Raphson method of root finding.

program main
    use legendre_stuff
    implicit none
    real (kind=8) :: x_i, root
    x_i = 0.1

    call newton_raphson(x_i, root, .false.)
    print *, 'root = ', root
    print *, 'Lambda(n, root) = ', Lambda(200, root)
end program main


subroutine newton_raphson(x_i, root, debug)
    use legendre_stuff
    implicit none
    integer :: iter
    logical, intent(in) :: debug
    real (kind=8), intent(in) :: x_i        ! initial guess value
    real (kind=8), intent(out) :: root
    real (kind=8) :: x_old, x_new
    real (kind=8) :: tol, rtol

    integer :: N = 200

    tol = 1.0d-16       ! absolute error tolerance
    rtol = 0.0d0        ! relative tolerance
    x_old = x_i

    if (debug) then
        print '(A3, 2A22, A25)', &
            'n', 'x_n', 'x_{n+1}', 'err = |x_{n+1} - x_n|'
    end if

    do iter = 1, 100
        x_new = x_old - Lambda(N, x_old) / Lambda_p(N, x_old)

        if (debug) then
            print '(I3, 2F22.16, E25.16)', &
                iter, x_old, x_new, abs(x_new - x_old)
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

