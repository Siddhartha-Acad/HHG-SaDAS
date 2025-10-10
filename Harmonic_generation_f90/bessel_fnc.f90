! Bessel function :: J_{n}(x)

! FORWARD RECURSION EQUATION
! J_{n}(x) = \frac{2(n-1)}{x} J_{n-1}(x) - J_{n-2}(x)

! BACKWARD RECURSION EQUATION
! J_{n-1}(x) = \frac{2n}{x} J_n(x) - J_{n+1}(x)


program main
    implicit none
    integer :: n_val, N = 27                 ! Forward recursion stable if (n < x)
    real (kind=8) :: x = 5.0d0
    real (kind=8) :: bessel_fwd,bessel_bwd
    external :: bessel_fwd, bessel_bwd

    open(unit=10, file='bessel_funcn.dat', status='replace', action='write')

    do n_val = 1, N
        write (10, *) n_val, bessel_fwd(n_val, x), &
            bessel_bwd(n_val, x), BESSEL_JN(n_val, x)
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
        Jk = (2.0d0 * real(k-1, kind=8) / x) * J1 - J0

        J0 = J1
        J1 = Jk
    end do

    bessel_fwd = Jk
end function bessel_fwd



pure real (kind=8) function bessel_bwd(n, x)
    implicit none
    integer :: k, Nmax
    integer, intent(in) :: n
    real (kind=8), intent(in) :: x
    real (kind=8) :: J_Nmax, J_Nmax_1      ! J_Nmax = J_{Nmax} ; J_Nmax_1 = J_{Nmax+1}
    real (kind=8) :: JN, Jk                ! JN = J_{N}(x) ; Jk = J_{N-1}(x)

    J_Nmax = 1.0d0
    J_Nmax_1 = 0.0d0

    Nmax = n + 100
    do k = Nmax, 1, -1
        Jk = (2.0d0 * real(k, kind=8) / x) * J_Nmax - J_Nmax_1

        if (k == n) then
            JN = J_Nmax
        end if

        J_Nmax_1 = J_Nmax
        J_Nmax = Jk
    end do

    JN = (BESSEL_J0(x) / Jk) * JN          ! norm_fact = BESSEL_J0(x) / Jk

    bessel_bwd = JN
end function bessel_bwd