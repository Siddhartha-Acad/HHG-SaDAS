! Newton-Raphson method of root finding.

program main
    implicit none
    real (kind=8) :: x_i
    real (kind=8) :: newton_raphson

    x_i = 3.0d0

    print *, newton_raphson(x_i)

end program main



pure function newton_raphson(x_i) result(root)
    implicit none
    real (kind=8), intent(in) :: x_i        ! initial guess value
    real (kind=8) :: root
    real (kind=8) :: x_n
    integer :: n

    do n = 1, 10
        x_n = x_i - sin(x_i) / cos(x_i)
        x_i = x_n
    end do

    root = x_i
end function newton_raphson