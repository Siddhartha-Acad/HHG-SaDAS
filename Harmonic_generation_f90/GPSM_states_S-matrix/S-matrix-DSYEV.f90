! fortran implementation for generating S-matrices.

program S_matrix_generator
    use parameters
    implicit none

#ifdef SAVE_STATES
    logical, parameter :: save_states = .true.
#else
    logical, parameter :: save_states = .false.
#endif

    integer :: i, j, l_val, recl_size
    real(kind=8), dimension(N-1) :: x, f_arr, fp_arr, d2_diag
    real(kind=8), dimension(N-1, N-1) :: fp_outer, d2_off_diag

    real(kind=8), dimension(N-1, N-1) :: H_matrix
    complex(kind=8), dimension(N-1, N-1) :: S_matrix

    character :: jobz, uplo
    integer :: lda, lwork, info
    real(kind=8), dimension(N-1) :: E_egval
    real(kind=8), allocatable :: work_array(:)

    character(len=8) :: l_str
    character(len=*), parameter :: green = char(27)//'[1;32m', &
                                   reset = char(27)//'[0m'


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


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                       Open relevant file(s)                        |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (.not. save_states) then
        inquire(iolength=recl_size) S_matrix
        print '(A, 1x, I0)', "Record length for S_matrix:", recl_size

        open(unit=20, file='data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
            form='unformatted', access='direct', recl=recl_size, status='replace')
    else
        inquire(iolength=recl_size) H_matrix(:, 1:check_n_states)
        print '(A, 1x, I0)', "Record length for eigenstates:", recl_size

        open(unit=30, file='data_GPSM_states_S-matrix/Eigenstates-DSYEV.bin', &
            form='unformatted', access='direct', recl=recl_size, status='replace')
    end if


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
            write(l_str, '(I0)') l_val
            if (.not. save_states) then
                S_matrix = matmul(H_matrix(:, 1:kmax) * &
                      spread(exp(cmplx(0.0d0, -E_egval(1:kmax) * dt / 2.0d0)), dim=1, ncopies=N-1), &
                      transpose(H_matrix(:, 1:kmax)))

                write(20, rec = l_val-m_qn+1) S_matrix                      ! rec = 1, 2, 3, ..., l_max
                print '(A, I2, A)', 'S(l=', l_val, ') matrix : ' // green // 'DONE' // reset

            else
                write(30, rec = l_val-m_qn+1) H_matrix(:, 1:check_n_states)
                print '(A, I2, A)', 'A(l=', l_val, ') states : ' // green // 'DONE' // reset
            end if
        else
            print '(A, I0)', 'DSYEV failed. info = ', info
        end if
    end do

    if (.not. save_states) close(20)
    if (save_states) close(30)
    deallocate(work_array)

end program S_matrix_generator

