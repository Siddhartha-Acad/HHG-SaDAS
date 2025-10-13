! compute all roots of a function in a given range.
! aiming to get rough estimation.

program main
    use legendre_stuff
    implicit none
    real (kind=8) :: dx
    integer :: i, root_count
    integer, parameter :: N = 5
    integer, parameter :: nop = 10000
    real (kind=8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real (kind=8), parameter :: xi = 0.0d0, xf = 1.0d0
    real (kind=8), allocatable :: roots(:)
    real (kind=8), dimension(nop) :: x, y

    dx = (xf - xi) / dble(nop - 1)

    do i = 1, nop
        x(i) = xi + (i - 1) * dx
    end do

    do i = 1, nop
        y(i) = -Lambda(N, x(i))**2
    end do

    root_count = 0
    allocate(roots(0))

    do i = 2, nop-1
        if (y(i-1) .lt. y(i) .and. y(i) .gt. y(i+1)) then
            root_count = root_count + 1
            roots = [roots, x(i)]
            print *, x(i)
        end if
    end do

    print *, allocated(roots)
    if (allocated(roots)) then
        deallocate(roots)
    end if

end program main
