! Newton-Raphson method of root finding.

program main
    use legendre_stuff
    implicit none
    real (kind=8) :: dx
    integer :: i, root_count
    integer, parameter :: N = 10
    integer, parameter :: nop = 1000
    real (kind=8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real (kind=8), parameter :: xi = 0.0d0, xf = 1.0d0
    real (kind=8), allocatable :: roots(:)
    real (kind=8), dimension(nop) :: x, y

    dx = (xf - xi) / dble(nop - 1)

    do i = 1, nop
        x(i) = xi + (i - 1) * dx
    end do

    do i = 1, nop
        y(i) = -Lambda(N, x(i))**2
    end do

    root_count = 0
    allocate(roots(0))

    do i = 2, nop-1
        if (y(i-1) .lt. y(i) .and. y(i) .gt. y(i+1)) then
            root_count = root_count + 1
            roots = [roots, x(i)]
        end if
    end do

    print '(A, I0, A)', '~~~~~~~~~~~~~~: Algo-3 :: N = ', N, ' :~~~~~~~~~~~~~~'
    if (mod(N, 2) .eq. 0 .and. size(roots) .eq. (N/2 - 1)) then
        print '(A, I2)', 'no. of initial guess values :', size(roots)
    else if (mod(N, 2) .ne. 0 .and. size(roots) .eq. (N-1)/2) then
        print '(A, I2)', 'no. of initial guess values :', size(roots)
    else
        print *, 'no. of initial guess values : wrong'
    end if

    if (allocated(roots)) then
        deallocate(roots)
    end if

end program main


subroutine newton_raphson(N, x_i, root, debug)
    use legendre_stuff
    implicit none
    integer :: iter
    integer, intent(in) :: N
    logical, intent(in) :: debug
    real (kind=8), intent(in)  :: x_i        ! initial guess value
    real (kind=8), intent(out) :: root
    real (kind=8) :: x_old, x_new
    real (kind=8) :: tol, rtol

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

