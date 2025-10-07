! legendre polynomial at a given x :: FORWARD RECURSION ALGORITHM

program main
    implicit none
    integer :: N
    real (kind=8) :: x(3) = [0.0d0, 0.5d0, 1.0d0]

    interface
        elemental real (kind=8) function legendre_p(n, x)
            integer, intent(in) :: n
            real (kind=8), intent(in) :: x
        end function legendre_p
    end interface

    N = 1
    print *, legendre_p(N, x)

end program main


elemental real (kind=8) function legendre_p(n, x)
    implicit none
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8) :: Pk
    real (kind=8) :: P0, P1
    integer :: k

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