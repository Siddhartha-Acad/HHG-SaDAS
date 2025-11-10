! File: Algo-3.f90
! Project: HHG-SaDAS
! 
! Author: Siddhartha Mithiya
! Affiliation: Indian Institute of Technology (IIT) Mandi
! License: MIT License
! Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git
! 
! Reference
! ------------
! Appendix A: "An efficient algorithm to numerically calculate the Gauss–Lobatto collocation points."
! 
! Notes
! ------------
! - Generates high-precision collocation points (accuracy <= O(10^-15))
! - This fortran code is written from scratch, no external dependencies.
! - This script is part of the HHG-SaDAS package: built for my Master of Science (Research) thesis:
!   "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
! --------------------------------------------------------------------------------

program main
    use timer_mod
    use legendre_stuff
    use parameters, only: green, yellow, reset
#ifndef N_VAL
    use parameters, only: N
#endif
    implicit none
#ifdef N_VAL
    integer, parameter :: N = N_VAL
    logical, parameter :: N_from_macro = .true.
#else
    logical, parameter :: N_from_macro = .false.
#endif

    integer :: i, root_count
    integer, parameter :: nop = 3000

    logical :: debug_newton = .false.
    logical :: print_colloc_pt = .false.
    real(kind=8), parameter :: xi_i = -1.0d0, xi_f = 1.0d0
    real(kind=8), dimension(N-1) :: colloc_pt
    real(kind=8), dimension(nop) :: x_map, y
    real(kind=8), allocatable :: roots(:)
    real(kind=8) :: xi, dx, guess
    real(kind=8) :: exec_time
    character(len=60) :: file_name

    call tick()                     ! start measuring time
    if (N_from_macro) then
        write(file_name, '(A, I0, A)') 'Algo-3_N=', N, '_Gauss_Lobatto_collocation_points.dat'
    else
        file_name = 'Algo-3_Gauss_Lobatto_collocation_points.dat'
    end if

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
        print '(A)', '  #        Initial Guess      Collocation point x(j)'
        print '(A)', ' ---  ---------------------  -----------------------'
    end if

    if (.not. ((mod(N, 2) .eq. 0 .and. root_count .eq. (N/2 - 1)) .or. &
               (mod(N, 2) .ne. 0 .and. root_count .eq. (N-1)/2))) then
        write(0, '(A,I0)') "ERROR: Incorrect number of initial guess values : increase 'nop', currently nop = ", nop
        stop 1
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

    print '(1x, A, I0)', 'N = ', N
    print '(1x, A, I0)', 'collocation point = ', N - 1

    if (mod(N, 2) .eq. 0) then
        colloc_pt(1:root_count) = -roots(root_count:1:-1)     ! First half: reversed negative roots
        colloc_pt(root_count + 1) = 0.0d0
        colloc_pt(root_count+2:N-1) = roots(1:root_count)     ! Second half: positive roots
    else
        colloc_pt(1:root_count) = -roots(root_count:1:-1)     ! First half: reversed negative roots
        colloc_pt(root_count+1:N-1) = roots(1:root_count)     ! Second half: positive roots
    end if

    if (allocated(roots)) then
        deallocate(roots)
    end if

    call tock(exec_time)                     ! stop measuring time

    open(unit=10, file=file_name, status='replace', action='write')
    write(10, *) colloc_pt
    close(10)

    print *
    print '(A, A, F0.5, A, A)', " Execution Wall-time: ", green, exec_time, reset, " sec"
    print '(A)', ' File created: ' // yellow // file_name // reset
    print *


contains
    pure real (kind=8) function f_rev(xi)
        implicit none
        real (kind=8), intent(in) :: xi
        real (kind=8), parameter :: L_map = 0.5d0, alpha = 1.0d0
        ! alpha = 2.0d0 * L_map

        f_rev = 1.0d0 - L_map * ((1-xi) / (1+xi+alpha))
    end function f_rev

end program main

