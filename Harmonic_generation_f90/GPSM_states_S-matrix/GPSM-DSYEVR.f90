! fortran implementation for solving 
! time independent Schrodinger equation using GPSM

program GPSM
    use parameters
    implicit none
    integer :: i, j
    real(kind=8), dimension(N-1) :: x                   ! collocation points
    real(kind=8), dimension(N-1, N-1) :: H_matrix
    
    character :: jobz, range, uplo
    real(kind=8) :: vl, vu, abstol
    real(kind=8), allocatable :: work_arr(:)
    real(kind=8), dimension(N-1) :: E_egval
    real(kind=8), dimension(N-1, kmax) :: E_vect
    integer, dimension(2*kmax) :: isuppz
    integer, allocatable :: iwork_arr(:)
    integer :: il, iu, lda, ldz, info, lwork, liwork, m
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                     reading collocation points                     |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    open(unit=10, file='./../collocation_points/generator/Algo-3_N=200_Gauss_Lobatto_collocation_points.dat', &
         status='old', action='read')
        read(10, *) x
    close(10)
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                     real symmetric [H] matrix                      |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do j = 1, N-1       ! Fill upper triangle (good cache access)
        do i = 1, j
            H_matrix(i, j) = H(l_qn, i, j)
        end do
    end do
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !          energy eigenvalues and eigenvectors : [H] matrix          |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! CALL DSYEVR(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z,
    !             LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO)
    
    lda = N-1        ! Leading dimension of H_matrix (number of rows in memory)
    ldz = N-1        ! Leading dimension of E_vect (number of rows in memory)
    jobz = 'V'       ! 'V' = compute eigenvectors; 'N' = eigenvalues only
    range = 'I'      ! 'I' = select eigenvalues by index (IL..IU)
    uplo = 'U'       ! 'U' = upper triangle of H_matrix is stored/used
    vl = 0.0d0       ! Lower bound of eigenvalues (used if RANGE='V', ignored here)
    vu = 0.0d0       ! Upper bound of eigenvalues (used if RANGE='V', ignored here)
    il = 1           ! Index of the smallest eigenvalue to compute (1 = first)
    iu = kmax        ! Index of the largest eigenvalue to compute (kmax = last wanted)
    abstol = 0.0d0   ! Absolute tolerance for convergence; 0 -> use machine precision
    
    ! Setting LWORK = -1 activates workspace query mode in LAPACK.
    lwork = -1
    liwork = -1
    allocate(work_arr(1), iwork_arr(1))
    call DSYEVR(jobz, range, uplo, N-1, H_matrix, lda, vl, vu, il, iu, abstol, &
                m, E_egval, E_vect, ldz, isuppz, work_arr, lwork, iwork_arr, liwork, info)
    lwork = int(work_arr(1))
    liwork = iwork_arr(1)
    deallocate(work_arr, iwork_arr)
    
    allocate(work_arr(lwork), iwork_arr(liwork))
    call DSYEVR(jobz, range, uplo, N-1, H_matrix, lda, vl, vu, il, iu, abstol, &
                m, E_egval, E_vect, ldz, isuppz, work_arr, lwork, iwork_arr, liwork, info)
    deallocate(work_arr, iwork_arr)
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !               Writing eigenvectors to an output file               |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (info .eq. 0) then
        do i = 1, kmax
            print '(A, I0, A, F20.16)', 'E(', i, ') =', E_egval(i)
        end do
        
        open(unit=11, file='data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
             form='unformatted', access='stream', status='replace')
        do i = 1, N-1
            write(11) f(x(i)), E_vect(i, :)
        end do
        close(11)
        print '(A)', 'GPSM Eigenvctors: data_GPSM_states_S-matrix/GPSM-DSYEVR_states.bin'
    else
        print '(A, I2)', 'DSYEVR failed. info = ', info
    end if
    
    
contains
    pure real(kind=8) function d2(i, j)
        integer, intent(in) :: i, j
        
        if (i .ne. j) then
            d2 = -2.0d0 / (x(i) - x(j))**2
        else
            d2 = -N*(N+1) / (3.0d0*(1.0d0 - x(i)**2))
        end if
    end function d2
    
    
    pure real(kind=8) function f(x_val)
        real(kind=8), intent(in) :: x_val
        
        f = Lmap * (1.0d0 + x_val) / (1.0d0 - x_val + alpha_map)
    end function f
    
    
    pure real(kind=8) function f_p(x_val)
        real(kind=8), intent(in) :: x_val
        
        f_p = Lmap * (alpha_map + 2.0d0) / (1.0d0 - x_val + alpha_map)**2
    end function f_p
    
    
    pure real(kind=8) function H(l_val, i, j)
        integer, intent(in) :: l_val, i, j
        real(kind=8) :: term1, term2
        
        term1 = -0.5d0 * (1.0d0 / f_p(x(i))) * d2(i, j) * (1.0d0 / f_p(x(j)))
        if (i .ne. j) then
            H = term1
        else
            term2 = l_val * (l_val + 1) / (2.0d0 * f(x(i)) ** 2) - 1.0d0 / f(x(i))
            H = term1 + term2
        end if
    end function H
    
end program GPSM

