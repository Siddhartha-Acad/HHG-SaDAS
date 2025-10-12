! compute all roots of a function in a given range.
! using Newton-Raphson method.

program main
    implicit none
    real (kind=8) :: dx
    integer :: i, root_count
    integer, parameter :: nop = 10000
    real (kind=8), parameter :: xi = 0.0d0, xf = 10.0d0
    real (kind=8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real (kind=8), dimension(nop) :: x, y
    real (kind=8), allocatable :: roots(:)
    dx = (xf - xi) / dble(nop - 1)

    x(1) = xi
    do i = 1, nop-1
        x(i+1) = x(i) + dx
    end do

    y = -sin(x)**2

    root_count = 0
    allocate(roots(0))

    do i = 2, nop-1
        if (y(i-1) .lt. y(i) .and. y(i) .gt. y(i+1)) then
            root_count = root_count + 1
            roots = [roots, x(i)]
            print *, x(i) / pi
        end if
    end do

    print *, allocated(roots)
    if (allocated(roots)) then
        deallocate(roots)
    end if

end program main