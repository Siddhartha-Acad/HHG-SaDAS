program check_Split_operator
    use parameters
    implicit none

    integer :: i, l_val, S_recl_size, A_recl_size
    complex(kind=8) :: S_matrix(N-1, N-1)

    real(kind=8) :: max_abs_err
    real(kind=8), dimension(N-1) :: r, A_tilde
    real(kind=8), dimension(N-1, check_n_states) :: A

    character(len=*), parameter :: green = char(27)//'[1;32m'  ! bold + green
    character(len=*), parameter :: red   = char(27)//'[1;31m'  ! bold + red
    character(len=*), parameter :: reset = char(27)//'[0m'     ! reset style

    inquire(iolength=S_recl_size) S_matrix
    inquire(iolength=A_recl_size) A

    open(unit=10, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/Eigenstates-DSYEV.bin', &
            form='unformatted', access='direct', recl=A_recl_size, status='old')
    open(unit=11, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
            form='unformatted', access='direct', recl=S_recl_size, status='old')

    do l_val = m_qn, l_max + m_qn

    read(10, rec=l_val - m_qn + 1) A
    read(11, rec=l_val - m_qn + 1) S_matrix        ! rec = 1, 2, 3, ..., l_max

    print '(A11, A, A11, A, A6)', 'S(l) @ A(n)', ' | ', 'abs_max_err', ' | ', 'status'
    print '(A)', repeat('-', 35)

    do i = 1, check_n_states
        A_tilde = matmul(S_matrix, cmplx(A(:, i), 0.0d0, kind=8))
        max_abs_err = maxval(abs(A_tilde - A(:, i)))

        if (max_abs_err .lt. 1.0d-3) then
            write(*,'(A,I0,A,I0,A,ES11.3,A,A)', advance='no') &
                'S(', m_qn, ') @ A(', i, ') | ', max_abs_err, ' | ', green//'[PASS]'//reset
            print *  ! newline after pass/fail
        else
            write(*,'(A,I0,A,I0,A,ES11.3,A,A)', advance='no') &
                'S(', m_qn, ') @ A(', i, ') | ', max_abs_err, ' | ', red//'[FAIL]'//reset
            print *  ! newline after pass/fail
        end if
    end do
    end do
    close(10)
    close(11)
end program check_Split_operator
