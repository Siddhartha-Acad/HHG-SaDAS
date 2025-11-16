! Gauss_Legendre.f90 :: golub–welsch algorithm

program main
    implicit none
    integer, parameter :: L = 5
    integer, parameter :: N = L + 1         ! to keep it simple.
    
    real(kind=8), dimension(N) :: diag, weights     ! diag will hold the nodes.
    real(kind=8), dimension(N-1) :: off_diag
    real(kind=8), dimension(N, N) :: Egvects
    
    integer :: i, info
    real(kind=8), allocatable :: work(:)
    
    diag = 0.0d0
    do i = 1, N-1
        off_diag(i) = dble(i) / sqrt(dble(4*i*i) - 1.0d0)
    end do

    allocate(work(max(1, 2*n-2)))
    
    ! SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
    call dstev('V', N, diag, off_diag, Egvects, N, work, info)

    if (info .ne. 0) then
     print *, "LAPACK dstev failed, INFO=", info
     stop
    end if

    ! LAPACK returns eigenvalues in diag(1:N) ascending, eigenvectors in columns of Egvects
    ! Weights: w_i = mu0 * (v_{0,i})^2. For Legendre on [-1,1], mu0 = 2 and
    ! v_{0,i} corresponds to Egvects(1,i) (first row, column i)
    do i = 1, N
     weights(i) = 2.0d0 * ( Egvects(1,i) * Egvects(1,i) )   ! = 2 * (v_{1,i})^2
    end do

    print *
    print '(2x,A,I0)', "Gauss-Legendre collocation points for L+1 = ", N
    print '(2A25)', "Node (x_i)", "Weight (w_i)"
    do i = 1, N
        print '(2(ES25.16))', diag(i), weights(i)
    end do
    
    deallocate(work)

end program main
