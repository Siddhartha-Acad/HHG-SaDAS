program read_GPSM_data
    use parameters
    implicit none
    integer :: i
    real :: rad_theory
    real(kind=8) :: rad_gpsm, rel_err
    real(kind=8), dimension(N-1) :: r
    real(kind=8), dimension(N-1, total_states) :: A

    open(unit=12, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
            form='unformatted', access='stream', status='old')
        read(12) r, A
    close(12)

    print '(A3, A, A11, A, A18, A, A10, A, A6)', ' n', ' | ', 'radius_Th', ' | ', 'radius_GPSM', ' | ', 'rel_err', ' | ', 'status'
    print '(A)', repeat('-', 62)

    do i = 1, total_states
        rad_theory = 0.5 * real(3*i**2 - l_qn*(l_qn+1))
        rad_gpsm = sum(r * A(:, i)**2)
        rel_err = abs(rad_gpsm - rad_theory) / rad_theory

        if (rel_err .lt. 1.0d-3) then
            print '(I3, A, F11.6, A, F18.13, A, ES10.2, A, A)', & 
                    i, ' | ', rad_theory, ' | ', rad_gpsm, ' | ', rel_err, ' | ', '[PASS]'
        else
            print '(I3, A, F11.6, A, F18.13, A, ES10.2, A, A)', &
                    i, ' | ', rad_theory, ' | ', rad_gpsm, ' | ', rel_err, ' | ', '[FAIL]'
        end if
    end do

end program read_GPSM_data
