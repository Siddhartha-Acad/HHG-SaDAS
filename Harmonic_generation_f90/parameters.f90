! parameters.f90

module parameters
    implicit none
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                GPSM Parameters                 |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer, parameter :: N = 200        ! P'_N(xj) = 0 ; radial grid size: len(colloc_pt) = N-1
    integer, parameter :: kmax = 5       ! number of GPSM states (maximum k index) in S matrix
    real(kind=8), parameter :: Lmap = 80.0d0, r_max = 200.0d0    ! radial mapping parameters
    real(kind=8), parameter :: alpha_map = 2.0d0 * Lmap / r_max


end module parameters
