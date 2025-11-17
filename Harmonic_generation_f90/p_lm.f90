
program test_legendre
    use legendre_stuff
    implicit none
    
    integer :: l, m
    real(8) :: x, result
    
    l = 3
    m = 2
    x = 0.5d0
    
    result = a_legendre(l, m, x)
    print *, "P_", l, "^", m, "(", x, ") = ", result
    
end program test_legendre
