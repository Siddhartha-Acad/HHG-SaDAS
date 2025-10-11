! compute all roots of a function in a given range.
! using Newton-Raphson method.

program main
    implicit none
    integer :: i
    real (kind=8) :: dx, max
    integer, parameter :: nop = 10000
    real (kind=8), parameter :: xi = 0.0d0, xf = 10.0d0
    real (kind=8), dimension(nop) :: x, y

    dx = (xf - xi) / dble(nop - 1)

    x(1) = xi
    do i = 1, nop-1
        x(i+1) = x(i) + dx
    end do

    y = sin(x)

    do i = 2, nop-1
        if (y(i-1) .lt. y(i) .and. y(i) .gt. y(i+1)) then
            print *, x(i) / 3.14159_8
        end if
    end do



end program main