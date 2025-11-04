! fortran implementation for generating S-matrices.

program S_matrix_generator
    use parameters
    implicit none
    integer :: i, j, l_val
    real(kind=8), dimension(N-1) :: x, f_arr, fp_arr, d2_diag
    real(kind=8), dimension(N-1, N-1) :: fp_outer, d2_off_diag

    real(kind=8), dimension(N-1, N-1) :: H_matrix
    complex(kind=8), dimension(N-1, N-1) :: S_matrix

    character :: jobz, uplo
    integer :: lda, lwork, info
    real(kind=8), dimension(N-1) :: E_egval
    real(kind=8), allocatable :: work_array(:)


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !      collocation points & Precompute l_val independent terms       |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    open(unit=10, file='./../collocation_points/generator/Algo-3_N=200_Gauss_Lobatto_collocation_points.dat', &
         status='old', action='read')
        read(10, *) x
    close(10)

    f_arr  = Lmap * (1.0d0 + x) / (1.0d0 - x + alpha_map)
    fp_arr = Lmap * (alpha_map + 2.0d0) / (1.0d0 - x + alpha_map)**2
    fp_outer = spread(fp_arr, dim=2, ncopies=N-1) * transpose(spread(fp_arr, dim=2, ncopies=N-1))

    d2_diag = -dble(N*(N+1)) / (3.0d0 * (1.0d0 - x**2))                    ! d(2)_ij : i = j
    d2_off_diag = 1.0d0 / (spread(x, dim=2, ncopies=N-1) - &
                           transpose(spread(x, dim=2, ncopies=N-1)))**2    ! d(2)_ij : i != j

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                   eigen-decomposition workspace                    |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
    lda = N-1        ! Leading dimension of H_matrix (number of rows in memory)
    jobz = 'V'       ! 'V' = compute eigenvectors; 'N' = eigenvalues only
    uplo = 'U'       ! 'U' = upper triangle of H_matrix is stored/used

    lwork = -1       ! Setting LWORK = -1 activates workspace query mode in LAPACK.
    allocate(work_array(1))
    call DSYEV(jobz, uplo, N-1, H_matrix, lda, E_egval, work_array, lwork, info)
    lwork = int(work_array(1)) 
    deallocate(work_array)

    allocate(work_array(lwork))
    open(unit=20, file='data_GPSM_states_S-matrix/S_matrices.bin', form='unformatted', access='sequential', status='replace')

    do l_val = m_qn, l_max + m_qn
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !                     real symmetric [H] matrix                      |
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        do j = 1, N-1
            ! Off-diagonal (Fill upper triangle)
            if (j > 1) then
                H_matrix(1:j-1, j) = d2_off_diag(1:j-1, j) / fp_outer(1:j-1, j)     ! -2 factor cancelled out.
            end if

            ! Diagonal element (i = j)
            H_matrix(j, j) = -0.5d0 * d2_diag(j) / fp_arr(j)**2 &
                            + dble(l_val*(l_val+1)) / (2.0d0 * f_arr(j)**2) - 1.0d0 / f_arr(j)
        end do

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !                      [H]-matrix & [S]-matrix                       |
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call DSYEV('V', 'U', N-1, H_matrix, N-1, E_egval, work_array, lwork, info)

        if (info .eq. 0) then
            S_matrix = matmul(H_matrix(:, 1:kmax) * &
                              spread(exp(cmplx(0.0d0, -E_egval(1:kmax) * dt / 2.0d0)), dim=1, ncopies=N-1), &
                              transpose(H_matrix(:, 1:kmax)))

            print '(A, I3, A)', 'S-matrix for l =', l_val, ' : DONE'
            write(20) S_matrix

        else
            print '(A, I2)', 'DSYEV failed. info = ', info
        end if
    end do

    close(20)

    deallocate(work_array)
end program S_matrix_generator

