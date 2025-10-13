! Newton-Raphson method of root finding.

program main
    use legendre_stuff
    implicit none
    integer :: i, root_count
    integer, parameter :: N = 10
    integer, parameter :: nop = 100
    real (kind=8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real (kind=8), parameter :: xi_i = -1.0d0, xi_f = 1.0d0
    real (kind=8), dimension(nop) :: x_map, y
    real (kind=8), dimension(N-1) :: colloc_pt
    real (kind=8), allocatable :: roots(:)
    real (kind=8), external :: f_rev
    real (kind=8) :: xi, dx, guess

    dx = (xi_f - xi_i) / dble(nop - 1)

    do i = 1, nop
        xi = xi_i + (i - 1) * dx
        x_map(i) = f_rev(xi)
        y(i) = -Lambda(N, x_map(i))**2
    end do

    root_count = 0
    allocate(roots(0))

    do i = 2, nop-1
        if (y(i-1) .lt. y(i) .and. y(i) .gt. y(i+1)) then
            root_count = root_count + 1
            roots = [roots, x_map(i)]
        end if
    end do

    print '(A)', ' '
    print '(A, I0, A)', '  ~~~~~~~~~~~~~~~: Algo-3 :: N = ', N, ' :~~~~~~~~~~~~~~~'
    print '(A)', ' '
    print '(A)', '   #        Initial Guess      Collocation point x(j)'
    print '(A)', '  ---  ---------------------  -----------------------'

    if (.not. ((mod(N, 2) .eq. 0 .and. root_count .eq. (N/2 - 1)) .or. &
               (mod(N, 2) .ne. 0 .and. root_count .eq. (N-1)/2))) then
        stop 'ERROR: Incorrect number of initial guess values'
    end if

    do i = 1, root_count
        guess = roots(i)
        call newton_raphson(N, roots(i), roots(i), .false.)
        print '(I4, 2X, F21.16, 2X, F23.16)', i, guess, roots(i)
    end do

    print '(A)', ' '
    print '(A, I0)', '  no. of collocation point : ', N - 1

    if (mod(N, 2) .eq. 0) then
        ! First half: reversed negative roots
        do i = 1, root_count
            colloc_pt(i) = -roots(root_count - i + 1)
        end do

        ! Middle element
        colloc_pt(root_count + 1) = 0.0d0

        ! Second half: positive roots
        do i = root_count + 2, N - 1
            colloc_pt(i) = roots(i - root_count - 1)
        end do
    else
        ! First half: reversed negative roots
        do i = 1, root_count
            colloc_pt(i) = -roots(root_count - i + 1)
        end do

        ! Second half: positive roots
        do i = root_count + 1, N - 1
            colloc_pt(i) = roots(i - root_count)
        end do
    end if

    if (allocated(roots)) then
        deallocate(roots)
    end if
end program main



pure real (kind=8) function f_rev(xi)
    implicit none
    real (kind=8), intent(in) :: xi
    real (kind=8) :: L_map, alpha

    L_map = 0.5d0
    alpha = 2.0d0 * L_map

    f_rev = 1.0d0 - L_map * ((1-xi) / (1+xi+alpha))
end function f_rev


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

