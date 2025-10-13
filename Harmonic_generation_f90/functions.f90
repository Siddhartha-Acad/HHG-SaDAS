! This script will contain the functions and subroutines
! that will be useful in other scripts...

module legendre_stuff
    implicit none

contains

pure real (kind=8) function legendre(n, x)
    ! FORWARD RECURSION ALGORITHM
    ! P_n(x) = \frac{1}{n}\left[ (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x) \right]

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


pure real (kind=8) function Lambda(n, x)
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x

    Lambda = legendre(n-1, x) - legendre(n+1, x)
end function Lambda


pure real (kind=8) function Lambda_p(n, x)
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8) :: coff_1, coff_2

    coff_1 = dble(n*(n-1)) / dble(2*n-1)
    coff_2 = dble((n+1)*(n+2)) / dble(2*n+3)

    Lambda_p = (coff_1*Lambda(n-1, x) - &
                coff_2*Lambda(n+1, x)) / (1.0d0 - x**2)
end function Lambda_p


end module legendre_stuff