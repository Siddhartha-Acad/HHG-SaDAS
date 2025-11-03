! parameters.f90

module parameters
    implicit none
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                      Atom                      |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer, parameter :: n_qn = 1, l_qn = 0, m_qn = 0        ! defines initial state. (qn = quantum number)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                GPSM Parameters                 |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer, parameter :: N = 200        ! P'_N(xj) = 0 ; radial grid size: len(colloc_pt) = N-1
    integer, parameter :: L = 20         ! must be >= l ; angular grid size: len(theta_k) = L+1
    integer, parameter :: l_max = 5     ! Number of partial waves = number of S-matrices = l_max+1
    integer, parameter :: kmax = 6       ! number of GPSM states (maximum k index) in S matrix
    real(kind=8), parameter :: Lmap = 80.0d0, r_max = 200.0d0    ! radial mapping parameters
    real(kind=8), parameter :: alpha_map = 2.0d0 * Lmap / r_max

    real(kind=8), parameter :: dt = 0.1


end module parameters
