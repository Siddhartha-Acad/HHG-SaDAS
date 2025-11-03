! fortran implementation for generating S-matrices.

program S_matrix_generator
    use parameters
    implicit none
    integer :: i, j, k, l_val
    real(kind=8), dimension(N-1) :: x                   ! collocation points
    real(kind=8), dimension(N-1, N-1) :: H_matrix
    complex(kind=8), dimension(N-1, N-1) :: S_matrix
    complex(kind=8) :: s_ij

    character :: jobz, uplo
    integer :: lda, lwork, info
    real(kind=8), dimension(N-1) :: E_egval
    real(kind=8), allocatable :: work_array(:)


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                     reading collocation points                     |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    open(unit=10, file='./../collocation_points/generator/Algo-3_N=200_Gauss_Lobatto_collocation_points.dat', &
         status='old', action='read')
        read(10, *) x
    close(10)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                   eigen-decomposition workspace                    |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
    lda = N-1        ! Leading dimension of H_matrix (number of rows in memory)
    jobz = 'V'       ! 'V' = compute eigenvectors; 'N' = eigenvalues only
    uplo = 'U'       ! 'U' = upper triangle of H_matrix is stored/used

    lwork = -1      ! Setting LWORK = -1 activates workspace query mode in LAPACK.
    allocate(work_array(1))
    call DSYEV(jobz, uplo, N-1, H_matrix, lda, E_egval, work_array, lwork, info)
    lwork = int(work_array(1)) 
    deallocate(work_array)

    allocate(work_array(lwork))

    do l_val = m_qn, l_max + m_qn
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !                     real symmetric [H] matrix                      |
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        do j = 1, N-1       ! Fill upper triangle (good cache access)
            do i = 1, j
                H_matrix(i, j) = H(l_val, i, j)
            end do
        end do

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !                      [H]-matrix & [S]-matrix                       |
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call DSYEV('V', 'U', N-1, H_matrix, N-1, E_egval, work_array, lwork, info)

        if (info .eq. 0) then
            S_matrix = matmul(H_matrix(:, 1:kmax) * &
                              spread(exp(cmplx(0.0d0, -E_egval(1:kmax) * dt / 2.0d0)), dim=1, ncopies=N-1), &
                              transpose(H_matrix(:, 1:kmax)))
            print '(I3, 5ES20.10)', l_val, real(S_matrix(1, 1:5))
        else
            print '(A, I2)', 'DSYEV failed. info = ', info
        end if
    end do

    deallocate(work_array)


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

end program S_matrix_generator

