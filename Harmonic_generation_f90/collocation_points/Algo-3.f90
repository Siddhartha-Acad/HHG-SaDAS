! Newton-Raphson method of root finding.

program main
    implicit none
    real (kind=8) :: x_i, root
    real (kind=8), external :: Lambda
    x_i = 0.1

    call newton_raphson(x_i, root, .false.)
    print *, 'root = ', root
    print *, 'Lambda(n, root) = ', Lambda(200, root)
end program main


pure real (kind=8) function legendre(n, x)
    ! FORWARD RECURSION ALGORITHM
    ! P_n(x) = \frac{1}{n}\left[ (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x) \right]

    implicit none
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8) :: P0, P1, Pk
    integer :: k

    P0 = 1.0d0
    P1 = x

    if (n==0) then
        legendre = 1.0d0
        return
    else if (n==1) then
        legendre = x
        return
    end if

    do k = 2, n
        Pk = ((2*k-1)*x*P1 - (k-1)*P0) / k
        P0 = P1
        P1 = Pk
    end do

    legendre = Pk
end function legendre


real (kind=8) function Lambda(n, x)
    implicit none
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8), external :: legendre

    Lambda = legendre(n-1, x) - legendre(n+1, x)
end function Lambda


real (kind=8) function Lambda_p(n, x)
    implicit none
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8), external :: Lambda
    real (kind=8) :: coff_1, coff_2

    coff_1 = dble(n*(n-1)) / dble(2*n-1)
    coff_2 = dble((n+1)*(n+2)) / dble(2*n+3)

    Lambda_p = (coff_1*Lambda(n-1, x) - &
                coff_2*Lambda(n+1, x)) / (1.0d0 - x**2)
end function Lambda_p


subroutine newton_raphson(x_i, root, debug)
    implicit none
    integer :: iter
    real (kind=8), intent(in) :: x_i        ! initial guess value
    real (kind=8), intent(out) :: root
    logical, intent(in) :: debug
    real (kind=8) :: x_old, x_new
    real (kind=8) :: tol, rtol
    real (kind=8), external :: Lambda, Lambda_p

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

