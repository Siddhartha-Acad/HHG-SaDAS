! Bessel function at a given x :: BACKWARD RECURSION ALGORITHM
! J_{n-1}(x) = \frac{2n}{x} J_n(x) - J_{n+1}(x)

program main
    implicit none
    integer :: k, N = 20
    integer :: Nmax
    real (kind=8) :: norm_fact
    real (kind=8) :: x = 50.0d0
    real (kind=8) :: J_Nmax, J_Nmax_1        ! J_Nmax = J_{Nmax} ; J_Nmax_1 = J_{Nmax+1}
    real (kind=8) :: JN, Jk, J0
    real (kind=8) :: JN_exact, rel_err

    J_Nmax = 1.0d0
    J_Nmax_1 = 0.0d0

    Nmax = N + 50

    do k = Nmax, 1, -1
        Jk = (2*k / x) * J_Nmax - J_Nmax_1

        if (k == N) then
            JN = J_Nmax
        end if

        J_Nmax_1 = J_Nmax
        J_Nmax = Jk

    end do

    norm_fact = BESSEL_J0(x) / Jk

    JN = norm_fact * JN
    JN_exact = BESSEL_JN(N, x)

    rel_err = abs(JN_exact - JN) / abs(JN_exact)

    print *, JN, JN_exact
    print *, int(-log10(rel_err))


end program main