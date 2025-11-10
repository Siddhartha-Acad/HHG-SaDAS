program read_GPSM_data
    use parameters
    implicit none
    integer :: i
    real :: rad_theory
    real(kind=8) :: rad_gpsm, rel_err
    real(kind=8), dimension(N-1) :: r
    real(kind=8), dimension(N-1, total_states) :: A

    character(len=*), parameter :: green = char(27)//'[1;32m', &
                                   red   = char(27)//'[1;31m', &
                                   reset = char(27)//'[0m'


    open(unit=12, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
            form='unformatted', access='stream', status='old')
        read(12) r, A
    close(12)

    print '(A3, A, A3, A, A11, A, A18, A, A10, A, A6)', &
            ' n', ' |', ' l', ' | ', 'radius_Th', ' | ', 'radius_GPSM', ' | ', 'rel_err', ' | ', 'status'
    print '(A)', repeat('-', 67)
    do i = 1, total_states
        rad_theory = 0.5 * real(3*(i + l_qn)**2 - l_qn*(l_qn+1))    ! n_qn = (l_qn + 1)
        rad_gpsm = sum(r * A(:, i)**2)
        rel_err = abs(rad_gpsm - rad_theory) / rad_theory

        if (rel_err .lt. 1.0d-3) then
            write(*,'(I3, A, I3, A, F11.6, A, F18.13, A, ES10.2, A, A)', advance='no') &
                (i + l_qn), ' |', l_qn, ' | ', rad_theory, ' | ', rad_gpsm, ' | ', rel_err, ' | ', green//'[PASS]'//reset
            print *   ! newline
        else
            write(*,'(I3, A, I3, A, F11.6, A, F18.13, A, ES10.2, A, A)', advance='no') &
                (i + l_qn), ' |', l_qn, ' | ', rad_theory, ' | ', rad_gpsm, ' | ', rel_err, ' | ', red//'[FAIL]'//reset
            print *   ! newline
        end if
    end do
    print *
    print '(A)', 'Note:'
    print '(A)', '  rel_err = |rad_gpsm - rad_theory| / rad_theory'
    print '(A)', '  ' // green // '[PASS]' // reset // ' : rel_err < 1.0d-3'
    print '(A)', '  ' // red   // '[FAIL]' // reset // ' : rel_err > 1.0d-3'
    print *

end program read_GPSM_data
