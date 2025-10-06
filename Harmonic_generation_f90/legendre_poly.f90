! computing legendre polynomial at a given x.

program legender_poly
    implicit none

    real (kind=8) :: x
    real (kind=8) :: Pk
    real (kind=8) :: P0, P1
    integer :: k, N = 10

    x = 0.5d0

    P0 = 1.0d0
    P1 = x

    do k = 2, N
        Pk = ((2*k-1)*x*P1 - (k-1)*P0) / k

        P0 = P1
        P1 = Pk
    end do

    print *, Pk

end program legender_poly
