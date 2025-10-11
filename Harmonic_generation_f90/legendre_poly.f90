! legendre polynomial at a given x :: FORWARD RECURSION ALGORITHM
! P_n(x) = \frac{1}{n}\left[ (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x) \right]

program main
    implicit none
    integer :: n_val, N = 500
    real (kind=8) :: x = 1.0d0
    real (kind=8), external :: legendre_fwd

    open(unit=10, file='legendre_poly.dat', status='replace', action='write')

    do n_val = 1, N
        write (10, *) n_val, legendre_fwd(n_val, x)
    end do

    close(10)
end program main


pure real (kind=8) function legendre_fwd(n, x)
    implicit none
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8) :: P0, P1, Pk
    integer :: k

    P0 = 1.0d0
    P1 = x

    if (n==0) then
        legendre_fwd = 1.0d0
        return
    else if (n==1) then
        legendre_fwd = x
        return
    end if

    do k = 2, n
        Pk = ((2*k-1)*x*P1 - (k-1)*P0) / k
        P0 = P1
        P1 = Pk
    end do

    legendre_fwd = Pk
end function legendre_fwd

