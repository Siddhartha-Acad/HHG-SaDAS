! vector_time_evolution.f90 : fortran implementation of ../Harmonic_generation_py/vector_time_evolution.py

program main
    use parameters
    implicit none

    integer :: S_recl_size
    real(kind=8), dimension(N) :: diag, weights
    real(kind=8), dimension(N-1) :: r
    real(kind=8), dimension(N-1, total_states) :: A

    inquire(iolength=S_recl_size) S_matrix

    open(unit=20, file="Gauss_Legendre_collocation_points_and_weights.dat", status="old", action="read")
        read(20, *) diag, weights
    close(20)

    open(unit=12, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin', &
            form='unformatted', access='stream', status='old')
        read(12) r, A
    close(12)

    open(unit=11, file='../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin', &
            form='unformatted', access='direct', recl=S_recl_size, status='old')
    close(11)


end program main
