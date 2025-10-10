! Newton-Raphson method of root finding.

program main
    implicit none
    real (kind=8) :: x_i
    real (kind=8) :: root, newton_raphson

    x_i = 3.0d0

    print '(a3, 3a20)', 'n', 'x_old', 'x_new', 'rel_err'
    root = newton_raphson(x_i)

end program main



function newton_raphson(x_i) result(root)
    implicit none
    real (kind=8), intent(in) :: x_i        ! initial guess value
    real (kind=8) :: root, x_old, x_new
    real (kind=8) :: tol, rel_err
    integer :: n

    tol = 1.0d-4
    rel_err = 1.0d0
    x_old = x_i

    do n = 1, 10
        x_new = x_old - sin(x_old) / cos(x_old)

        rel_err = abs(x_new - x_old) / abs(x_new)

        print '(I3, 2F20.16, E25.16)', n, x_old, x_new, rel_err

        if (rel_err < tol) then
            root = x_new
            return
        end if

        x_old = x_new
    end do

    root = x_old
end function newton_raphson