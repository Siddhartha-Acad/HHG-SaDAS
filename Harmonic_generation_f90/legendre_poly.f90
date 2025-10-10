! legendre polynomial at a given x :: FORWARD RECURSION ALGORITHM
! P_n(x) = \frac{1}{n}\left[ (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x) \right]

program main
    implicit none
    integer :: n
    real (kind=8) :: x = 0.5d0
    real (kind=8), external :: legendre_fwd, legendre_bwd

    open(unit=10, file='legendre_poly.dat', status='replace', action='write')

    do n = 1, 500
        write (10, *) n, legendre_bwd(n, x)
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


pure real (kind=8) function legendre_bwd(n, x)
    implicit none
    integer :: k, Nmax
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8) :: P_Nmax, P_Nmax_1      ! P_Nmax = P_{Nmax} ; P_Nmax_1 = P_{Nmax+1}
    real (kind=8) :: PN, Pk                ! PN = P_{N}(x) ; Pk = P_{N-1}(x)

    P_Nmax = 1.0d0
    P_Nmax_1 = 0.0d0

    Nmax = n + 100
    do k = Nmax, 1, -1
        Pk = ((2.0d0 * real(k, kind=8) + 1.0d0) * x * P_Nmax - &
             (real(k, kind=8)+1.0d0) * P_Nmax_1) / real(k, kind=8)

        if (k == n) then
            PN = P_Nmax
        end if

        P_Nmax_1 = P_Nmax
        P_Nmax = Pk
    end do

    PN = (1.0d0 / Pk) * PN          ! norm_fact = P_0(x) / Pk

    legendre_bwd = PN
end function legendre_bwd