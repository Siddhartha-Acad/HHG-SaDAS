! File: Algo-3.f90
! Project: HHG-SaDAS
!
! $ cd ./Harmonic_generation_f90/collocation_points/
! $ gfortran -J.. -c ../functions.f90
! $ gfortran -I.. -c ./newton_raphson.f90
! $ gfortran -I.. ./Algo-3.f90 ./newton_raphson.o ./functions.o -o Algo-3.exe
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
    use parameters
    use legendre_stuff
    implicit none
    integer :: i, root_count
    integer, parameter :: nop = 3000

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
