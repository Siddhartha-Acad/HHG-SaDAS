! legendre polynomial at a given x.

program legender_poly
    implicit none
    integer :: N
    real (kind=8) :: x
    real (kind=8) :: legendre_p

    N = 0
    x = 0.5d0

    print *, legendre_p(N, x)

end program legender_poly


real (kind=8) function legendre_p(n, x)
    implicit none
    real (kind=8) :: x
    real (kind=8) :: Pk
    real (kind=8) :: P0, P1
    integer :: k, n

    P0 = 1.0d0
    P1 = x

    if (n==0) then
        legendre_p = 1.0d0
        return
    else if (n==1) then
        legendre_p = x
        return
    end if

    do k = 2, n
        Pk = ((2*k-1)*x*P1 - (k-1)*P0) / k
        P0 = P1
        P1 = Pk
    end do

    legendre_p = Pk
end function legendre_p