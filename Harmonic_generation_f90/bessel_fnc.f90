! Bessel function at a given x :: FORWARD RECURSION ALGORITHM

program main
    implicit none
    integer :: n_val, N = 20     ! Forward recursion stable if (n < x)
    real (kind=8) :: x = 1.0d0
    real (kind=8) :: bessel_fwd
    external :: bessel_fwd

    open(unit=10, file='bessel_fwd.dat', status='replace', action='write')

    do n_val = 1, N
        write (10, *) n_val, bessel_fwd(n_val, x), BESSEL_JN(n_val, x)
    end do

    close(10)
end program main


pure real (kind=8) function bessel_fwd(n, x)
    implicit none
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8) :: J0, J1, Jk
    integer :: k

    J0 = BESSEL_J0(x)
    J1 = BESSEL_J1(x)

    if (n==0) then
        bessel_fwd = BESSEL_J0(x)
        return
    else if (n==1) then
        bessel_fwd = BESSEL_J1(x)
        return
    end if

    do k = 2, n
        Jk = (2*(k-1) / x) * J1 - J0

        J0 = J1
        J1 = Jk
    end do

    bessel_fwd = Jk
end function bessel_fwd