! Subroutine for computing roots newton raphson method.

subroutine newton_raphson(N, x_i, root, debug)
    use legendre_stuff
    implicit none
    integer, intent(in) :: N
    logical, intent(in) :: debug
    real (kind=8), intent(in)  :: x_i        ! initial guess value
    real (kind=8), intent(out) :: root

    integer :: iter
    real (kind=8) :: x_old, x_new
    real (kind=8), parameter :: tol = 1.0d-15, rtol = 0.0d0

    ! tol = absolute tolerance
    ! rtol = relative tolerance

    if (debug) then
        print '(A)'
        print '(A)', '+-----+----------------------+----------------------+---------------------------+'
        print '(A)', '|  n  |         x_n          |       x_{n+1}        |   err = |x_{n+1} - x_n|   |'
        print '(A)', '+-----+----------------------+----------------------+---------------------------+'
    end if

    x_old = x_i

    do iter = 1, 50
        x_new = x_old - Lambda(N, x_old) / Lambda_p(N, x_old)

        if (debug) then
            print '(A, I3, A, F20.16, A, F20.16, A, E24.16, A)', &
                    '| ', iter, ' | ', x_old, ' | ', x_new, ' | ', abs(x_new - x_old), '  |'
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
