! compute all roots of a function in a given range.
! aiming to get rough estimation.

program main
    implicit none
    real (kind=8) :: dx
    integer :: i, root_count
    integer, parameter :: nop = 10000
    real (kind=8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real (kind=8), parameter :: xi = 0.0d0, xf = 1.0d0
    real (kind=8), allocatable :: roots(:)
    real (kind=8), dimension(nop) :: x, y

    integer, parameter :: N = 5
    real (kind=8), external :: Lambda

    dx = (xf - xi) / dble(nop - 1)
    x = [(xi + (i-1)*dx, i=1,nop)]

    do i = 1, nop
        y(i) = -Lambda(N, x(i))**2
    end do

    root_count = 0
    allocate(roots(0))

    do i = 2, nop-1
        if (y(i-1) .lt. y(i) .and. y(i) .gt. y(i+1)) then
            root_count = root_count + 1
            roots = [roots, x(i)]
            print *, x(i)
        end if
    end do

    print *, allocated(roots)
    if (allocated(roots)) then
        deallocate(roots)
    end if

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
        Pk = ((2*k-1.0d0)*x*P1 - (k-1.0d0)*P0) / k
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