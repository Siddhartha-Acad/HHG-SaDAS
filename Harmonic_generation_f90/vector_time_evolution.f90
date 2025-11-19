! vector_time_evolution.f90 : fortran implementation of ../Harmonic_generation_py/vector_time_evolution.py

program main
    use parameters
    implicit none

    integer :: S_recl_size
    integer(kind=8) :: A_pos_offset
    real(kind=8), dimension(N) :: roots, weights
    real(kind=8), dimension(N-1) :: r, A_r

    inquire(iolength=S_recl_size) S_matrix

    open(unit=20, file="Gauss_Legendre_collocation_points_and_weights.dat", status="old", action="read")
        read(20, *) roots, weights
    close(20)

    open(unit=12, file='GPSM-DSYEV_states.bin', access='stream', form='unformatted')
        read(12) r
        A_pos_offset = 8_int64 * ( (N-1) + (n_qn-1)*(N-1) )       ! Compute byte offset for column n
        read(12, pos=A_pos_offset) A_col                          ! Jump directly to column n
    close(12)

    open(unit=11, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
            form='unformatted', access='direct', recl=S_recl_size, status='old')
    close(11)


end program main
