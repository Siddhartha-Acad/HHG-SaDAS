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
    integer, parameter :: N = 200            ! P'_N(xj) = 0 ; radial grid size: len(colloc_pt) = N-1
    integer, parameter :: L = 20             ! must be >= l ; angular grid size: len(theta_k) = L+1
    integer, parameter :: l_max = 20         ! Number of partial waves = number of S-matrices = l_max+1
    integer, parameter :: kmax = 50          ! number of GPSM states (maximum k index) in S matrix
    integer, parameter :: total_states = 5   ! how many states you want to keep in the GPSM_state file (.dat)
    real(kind=8), parameter :: Lmap = 80.0d0, r_max = 200.0d0, r0 = 150.0d0    ! radial mapping parameters and absorber radius
    real(kind=8), parameter :: alpha_map = 2.0d0 * Lmap / r_max


    real(kind=8), parameter :: dt = 0.1
    integer, parameter :: check_n_states = 6    ! how many states you want to include in check_Split_operator.f90?

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                  ANSI colors                   |
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    character(len=*), parameter :: green  = char(27)//'[1;32m', &
                                   red    = char(27)//'[1;31m', &
                                   yellow = char(27)//'[1;33m', &
                                   white  = char(27)//'[1;37m', &
                                   reset  = char(27)//'[0m'

end module parameters
