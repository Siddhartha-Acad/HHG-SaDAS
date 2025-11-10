program check_Split_operator
    use timer_mod
    use parameters
    implicit none

    integer :: i, l_val, S_recl_size, A_recl_size
    complex(kind=8) :: S_matrix(N-1, N-1)

    real(kind=8) :: max_abs_err
    real(kind=8), dimension(N-1) :: A_tilde
    real(kind=8), dimension(N-1, check_n_states) :: A
    real(kind=8) :: exec_time

    integer, parameter :: CW = 8         ! column width (adjust if you want wider)
    character(len=CW) :: col             ! temporary column field
    character(len=8)  :: inum            ! enough for the integer -> text
    character(len=20) :: fmtA            ! format string like '(A5)'

    call tick()             ! start measuring time

    inquire(iolength=S_recl_size) S_matrix
    inquire(iolength=A_recl_size) A

    open(unit=10, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/Eigenstates-DSYEV.bin', &
            form='unformatted', access='direct', recl=A_recl_size, status='old')
    open(unit=11, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
            form='unformatted', access='direct', recl=S_recl_size, status='old')


    write(fmtA, '(A,I0,A)') '(A', CW, ')'                       ! format string dynamically based on CW
    print '(1x, A)', repeat('-', (1 + check_n_states) * (CW + 1))   ! Horizontal rule

    write(*, fmtA, advance='no') 'S @ A'
    do i = 1, check_n_states
        write(inum,'(I0)') i
        write(*,'(1X'//trim(fmtA)//')', advance='no') 'A(' // trim(inum) // ')'
    end do
    print *
    print '(1x, A)', repeat('-', (1 + check_n_states) * (CW + 1))


    do l_val = m_qn, l_max + m_qn
        read(10, rec=l_val - m_qn + 1) A
        read(11, rec=l_val - m_qn + 1) S_matrix

        write(inum,'(I0)') l_val
        write(*, '('//trim(fmtA)//',1X)', advance='no') 'S(' // trim(inum) // ')'

        do i = 1, check_n_states
            A_tilde = matmul(S_matrix, cmplx(A(:, i), 0.0d0, kind=8))
            max_abs_err = maxval(abs(A_tilde - A(:, i)))
            if (max_abs_err < 1.0d-3) then
                write(col, fmtA) '[PASS]'
                write(*,'(1X,A)', advance='no') green // trim(col) // reset
            else
                write(col, fmtA) '[FAIL]'
                write(*,'(1X,A)', advance='no') red // trim(col) // reset
            end if
        end do
        print *
    end do

    close(10)
    close(11)

    call tock(exec_time)            ! stop measuring time

    print '(1x, A)', repeat('-', (1 + check_n_states) * (CW + 1))
    print *
    print '(1x, A, A, F0.5, A, A)', "Execution Wall-time: ", green, exec_time, reset, " sec"
    print '(1x, A)', 'Note:'
    print '(1x, A)', '  err = max( | S(l) @ A(n, l) - A(n, l) | )'
    print '(1x, A)', '  ' // green // '[PASS]' // reset // ' : err < 1.0d-3'
    print '(1x, A)', '  ' // red   // '[FAIL]' // reset // ' : err > 1.0d-3'
    print '(1x, A)', repeat('-', (1 + check_n_states) * (CW + 1))

end program check_Split_operator
