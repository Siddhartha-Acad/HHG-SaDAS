! File: Algo-3.f90
! Project: HHG-SaDAS
! Code Description:
!     ! *** [Main Gauss-Lobatto collocation point generating code] ***
!     !
!     ! Following Appendix-A of my thesis:
!     ! - Algo-3 uses the non-equispaced grid x⁺(ξ) ∈ (0, 1), as in Eq.A.11,
!     !   to determine the roots of Λ_N(x) from the local maxima of -Λ_N(x)² (serving as the initial guesses).
!     ! - It calculates only half of the roots (those in the positive interval), while the other half
!     !   (in the negative interval) are obtained using the parity relation in Eq.A.10.
!     ! - This algorithm is graphically presented in the flowchart of Fig.A.4.
!
! >>> cd .\Harmonic_generation_f90\collocation_points\
! >>> gfortran -J.. -c ..\functions.f90
! >>> gfortran -I.. .\Algo-3.f90 .\functions.o -o Algo-3.exe
!
! Author: Siddhartha Mithiya
! Affiliation: Indian Institute of Technology (IIT) Mandi
! License: MIT License
! Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git
! 
! --------------------------------------------------------------------------------
! Notes:
! - Generates high-precision collocation points (accuracy <= O(10^-15))
! - This fortran code is written from scratch, no external dependencies.
! - This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
!   "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
! --------------------------------------------------------------------------------

program main
    use legendre_stuff
    implicit none
    integer :: i, root_count
    integer, parameter :: N = 10
    integer, parameter :: nop = 80

    logical :: debug_newton = .false.
    logical :: print_colloc_pt = .false.
    real (kind=8), parameter :: xi_i = -1.0d0, xi_f = 1.0d0
    real (kind=8), dimension(N-1) :: colloc_pt
    real (kind=8), dimension(nop) :: x_map, y
    real (kind=8), allocatable :: roots(:)
    real (kind=8), external :: f_rev
    real (kind=8) :: xi, dx, guess
    character(len=60) :: file_name

    dx = (xi_f - xi_i) / dble(nop - 1)

    do i = 1, nop
        xi = xi_i + (i - 1) * dx
        x_map(i) = f_rev(xi)
        y(i) = -Lambda(N, x_map(i))**2
    end do

    root_count = 0
    allocate(roots(N/2))    ! optimum length to hold roots.
                            ! N/2 > (N-1)/2 > N/2-1

    do i = 2, nop-1
        if (y(i-1) .lt. y(i) .and. y(i) .gt. y(i+1)) then
            root_count = root_count + 1
            roots(root_count) = x_map(i)
        end if
    end do

    roots = roots(:root_count)

    if (debug_newton .eqv. print_colloc_pt) then
        print_colloc_pt = .false.
    end if

    if (print_colloc_pt) then
        print '(A)', ' '
        print '(A)', '   #        Initial Guess      Collocation point x(j)'
        print '(A)', '  ---  ---------------------  -----------------------'
    end if

    if (.not. ((mod(N, 2) .eq. 0 .and. root_count .eq. (N/2 - 1)) .or. &
               (mod(N, 2) .ne. 0 .and. root_count .eq. (N-1)/2))) then
        stop 'ERROR: Incorrect number of initial guess values'
    end if

    do i = 1, root_count
        if (print_colloc_pt) then
            guess = roots(i)
        end if
        call newton_raphson(N, roots(i), roots(i), debug_newton)
        if (print_colloc_pt) then
            print '(I4, 2X, F21.16, 2X, F23.16)', i, guess, roots(i)
        end if
    end do

    print '(A)', ' '
    print '(A, I0)', '  No. of collocation point: ', N - 1

    if (mod(N, 2) .eq. 0) then
        ! First half: reversed negative roots
        do i = 1, root_count
            colloc_pt(i) = -roots(root_count - i + 1)
        end do

        colloc_pt(root_count + 1) = 0.0d0

        ! Second half: positive roots
        do i = root_count + 2, N - 1
            colloc_pt(i) = roots(i - root_count - 1)
        end do
    else
        ! First half: reversed negative roots
        do i = 1, root_count
            colloc_pt(i) = -roots(root_count - i + 1)
        end do

        ! Second half: positive roots
        do i = root_count + 1, N - 1
            colloc_pt(i) = roots(i - root_count)
        end do
    end if

    if (allocated(roots)) then
        deallocate(roots)
    end if

    write(file_name, '(A, I0, A)') 'Algo-3_N=', N, '_Gauss_Lobatto_collocation_points.dat'
    open(unit=10, file=file_name, status='replace', action='write')

    do i = 1, N-1
        write(10, *) colloc_pt(i)
    end do

    close(10)
    print '(2A)', '  File created: ', file_name; print '(A)', ' '
end program main



pure real (kind=8) function f_rev(xi)
    implicit none
    real (kind=8), intent(in) :: xi
    real (kind=8), parameter :: L_map = 0.5d0, alpha = 1.0d0
    ! alpha = 2.0d0 * L_map

    f_rev = 1.0d0 - L_map * ((1-xi) / (1+xi+alpha))
end function f_rev


subroutine newton_raphson(N, x_i, root, debug)
    use legendre_stuff
    implicit none
    integer, intent(in) :: N
    logical, intent(in) :: debug
    real (kind=8), intent(in)  :: x_i        ! initial guess value
    real (kind=8), intent(out) :: root

    integer :: iter
    real (kind=8) :: x_old, x_new
    real (kind=8), parameter :: tol = 1.0d-15, rtol = 0.0d0

    ! tol = absolute error tolerance
    ! rtol = relative tolerance

    if (debug) then
        print '(A)'
        print '(A)', '+-----+----------------------+----------------------+---------------------------+'
        print '(A)', '|  n  |         x_n          |       x_{n+1}        |   err = |x_{n+1} - x_n|   |'
        print '(A)', '+-----+----------------------+----------------------+---------------------------+'
    end if

    x_old = x_i

    do iter = 1, 50
        x_new = x_old - Lambda(N, x_old) / Lambda_p(N, x_old)

        if (debug) then
            print '(A, I3, A, F20.16, A, F20.16, A, E24.16, A)', &
                    '| ', iter, ' | ', x_old, ' | ', x_new, ' | ', abs(x_new - x_old), '  |'
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
