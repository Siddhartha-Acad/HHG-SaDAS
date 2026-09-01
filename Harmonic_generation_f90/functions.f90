! This script will contain the functions and subroutines
! that will be useful in other scripts...

module legendre_stuff
    implicit none

contains

    pure real(kind=8) function legendre(n, x)
        ! FORWARD RECURSION : Bonnet's recursion formula
        ! P_n(x) = \frac{1}{n}\left[ (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x) \right]
        
        integer :: k
        integer, intent(in) :: n
        real(kind=8), intent(in) :: x
        real(kind=8) :: P0, P1, Pk
        
        if (n==0) then
            legendre = 1.0d0
            return
        else if (n==1) then
            legendre = x
            return
        end if
        
        P0 = 1.0d0
        P1 = x
        
        do k = 2, n
            Pk = ((2*k-1)*x*P1 - (k-1)*P0) / dble(k)
            P0 = P1
            P1 = Pk
        end do
        
        legendre = Pk
    end function legendre


    pure real(kind=8) function Lambda(n, x)
        integer, intent(in) :: n
        real(kind=8), intent(in) :: x
        
        Lambda = legendre(n-1, x) - legendre(n+1, x)
    end function Lambda


    pure real(kind=8) function Lambda_p(n, x)
        integer, intent(in) :: n
        real(kind=8), intent(in) :: x
        real(kind=8) :: coff_1, coff_2
        
        coff_1 = dble(n*(n-1)) / dble(2*n-1)
        coff_2 = dble((n+1)*(n+2)) / dble(2*n+3)
        
        Lambda_p = (coff_1*Lambda(n-1, x) - &
                    coff_2*Lambda(n+1, x)) / (1.0d0 - x**2)
    end function Lambda_p


    real(kind=8) function a_legendre(l, m, x)
    ! associated Legendre polynomial P_l^m(x)
    ! Input: l (degree), m (order), x (argument, -1 <= x <= 1)
    ! Output: P_l^m(x)
        integer, intent(in) :: l, m
        real(kind=8), intent(in) :: x
        real(kind=8) :: pmm, pmmp1, pll
        real(kind=8) :: fact
        integer :: i
        
        ! Check validity
        if (m .lt. 0 .or. m .gt. l .or. abs(x) .gt. 1.0d0) then
            print *, "Error: Invalid parameters"
            a_legendre = 0.0d0
            return
        end if
        
        ! Diagonal Terms (no CS phase) :: P_m^m(x) = (2m-1)!! (1-x^2)^{m/2}
        pmm = 1.0d0
        if (m .gt. 0) then
            fact = 1.0d0
            do i = 1, m
                pmm = pmm * fact * sqrt(1.0d0 - x**2)
                fact = fact + 2.0d0
            end do
        end if
        
        if (l .eq. m) then
            a_legendre = pmm
            return
        end if
        
        ! One Step Off-Diagonal :: P_{m+1}^m(x) = x(2m+1) P_m^m(x)
        pmmp1 = x * (2*m + 1) * pmm
        if (l .eq. m+1) then
            a_legendre = pmmp1
            return
        end if
        
        ! General Three-Term Recurrence for l > m+1 :: P_l^m(x) = \frac{x(2l-1) P_{l-1}^m(x) - (l+m-1) P_{l-2}^m(x)}{l-m}
        do i = m + 2, l
            pll = (x * (2*i - 1) * pmmp1 - (i+ m - 1) * pmm) / (i - m)
            pmm = pmmp1
            pmmp1 = pll
        end do
        
        a_legendre = pmmp1
    end function a_legendre

end module legendre_stuff



subroutine newton_raphson(N, x_i, root, debug)
    use legendre_stuff
    implicit none
    integer, intent(in) :: N
    logical, intent(in) :: debug
    real(kind=8), intent(in)  :: x_i        ! initial guess value
    real(kind=8), intent(out) :: root
    
    integer :: iter
    real(kind=8) :: x_old, x_new
    real(kind=8), parameter :: tol = 1.0d-15, rtol = 0.0d0
    
    ! tol = absolute tolerance
    ! rtol = relative tolerance
    
    if (debug) then
        print '(A)'
        print '(1x, A)', '+-----+----------------------+----------------------+--------------------------+'
        print '(1x, A)', '|  n  |         x_n          |       x_{n+1}        |   err = |x_{n+1} - x_n|  |'
        print '(1x, A)', '+-----+----------------------+----------------------+--------------------------+'
    end if
    
    x_old = x_i
    do iter = 1, 50
        x_new = x_old - Lambda(N, x_old) / Lambda_p(N, x_old)
        
        if (debug) then
            print '(1x, A, I3, A, F20.16, A, F20.16, A, E24.16, A)', &
                    '| ', iter, ' | ', x_old, ' | ', x_new, ' | ', abs(x_new - x_old), ' |'
        end if
        
        if (abs(x_new - x_old) .lt. (tol + rtol*abs(x_new))) then
            root = x_new
            return
        end if
        
        x_old = x_new
    end do
    
    ! Reaching here is signature that convergence failed.
    root = 1.0d2
end subroutine newton_raphson


integer(kind=8) function factorial(n)
    integer, intent(in) :: n
    integer(kind=8) :: f
    integer :: i

    f = 1_8
    do i = 2, n
        f = f * int(i, kind=8)
    end do
    factorial = f
end function factorial


real(kind=8) function N_fact(l_val, m_val)
    implicit none
    integer, intent(in) :: l_val, m_val
    integer(kind=8), external :: factorial
    real(kind=8), parameter :: pi = acos(-1.0d0)
    real(kind=8) :: num, den
    integer :: CS_phase

    if (mod(abs(m_val), 2) == 0) then
        CS_phase = 1
    else
        CS_phase = -1
    end if

    num = dble(factorial(l_val - m_val))
    den = dble(factorial(l_val + m_val))

    N_fact = dble(CS_phase) * sqrt((2.0d0 * dble(l_val) + 1.0d0) * num / (4.0d0 * pi * den))
end function N_fact

real(kind=8) function C_fact(l_val, m_val)
    implicit none
    integer, intent(in) :: l_val, m_val
    integer(kind=8), external :: factorial
    if (l_val < 0 .or. m_val < 0 .or. m_val > l_val) then
        C_fact = 0.0d0
        return
    end if
    C_fact = 2.0d0 * dble(factorial(l_val + m_val)) / ((2.0d0 * dble(l_val) + 1.0d0) * dble(factorial(l_val - m_val)))
end function C_fact

real(kind=8) function Y_lm(l_val, m_val, x)
    use legendre_stuff, only: a_legendre
    implicit none
    integer, intent(in) :: l_val, m_val
    real(kind=8), intent(in) :: x
    real(kind=8), external :: N_fact

    Y_lm = N_fact(l_val, m_val) * a_legendre(l_val, m_val, x)
end function Y_lm

