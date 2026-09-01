! vector_time_evolution.f90 : Fortran implementation of
! ../Harmonic_generation_py/vector_time_evolution.py
!
! This version reproduces the short validation loop used in the Python
! reference: apply S(l), reconstruct psi(theta,r), multiply by the interaction
! phase, project back to g_lm(r), and absorb the outer region.
!
! --------------------------------------------------------------------------
! [PERFORMANCE FIXES vs. the original version]
!   1) a_legendre(l,0,roots(j)) was being re-evaluated inside the ti-loop,
!      (N-1)*l_max*(L+1) times PER STEP, even though the exact same values
!      were already computed once into Y_T before the loop. Now both Y_T
!      and P_table are built ONCE before the time loop from a single
!      pass over a_legendre, and the loop only ever reads P_table.
!   2) norm_factor = N_fact(l-1,0)*C_fact(l-1,0) was recomputed N-1 times
!      per l per step even though it does not depend on i. Now precomputed
!      once into norm_factor(l_max) before the loop.
!   3) S_matrix_all, init_glm and glm_tilde were declared with the l-index
!      as the LEADING (fastest-varying) Fortran dimension, so slices like
!      S_matrix_all(l_idx,:,:) and init_glm(l_idx,:) used in matmul() were
!      non-contiguous. gfortran must copy non-contiguous arguments into a
!      temporary before calling the (intrinsic) matmul, on every call. These
!      three arrays are now dimensioned with l as the TRAILING dimension so
!      S_matrix_all(:,:,l_idx) / init_glm(:,l_idx) / glm_tilde(:,l_idx) are
!      contiguous columns/blocks and matmul() gets them with no hidden copy.
!   4) Y_T is stored with l as the leading dimension (matches how it is
!      walked in the innermost loop of the angular reconstruction), and
!      P_table / weighted_P are stored with j leading (matches how they are
!      walked in the innermost loop of the projection step) -- both for
!      cache-friendly, unit-stride access in the manual reduction loops.
!   5) OpenMP added on top of the above. All four per-step hot spots inside
!      the ti-loop are embarrassingly parallel over their OUTER index
!      (l_idx for steps 1/3/4, j/i for step 2) because each outer-index
!      iteration reads its own slice of the inputs and writes to a disjoint
!      slice of the outputs -- no cross-iteration dependency, so no reduction
!      clause or critical section is needed anywhere.
!
!   [FURTHER OPTIMIZATIONS -- merged in from the "chatgpt" version, aimed
!    purely at making the already-correct algorithm above run faster]
!
!   6) ONE PERSISTENT OpenMP PARALLEL REGION around the *entire* time loop.
!      The previous version opened and closed a `!$omp parallel do` team
!      FOUR separate times per timestep (steps 1-4), which means the OpenMP
!      runtime forked/joined the thread team 4*time_step times over a run.
!      Now a single `!$omp parallel` region is opened once, before `do ti =
!      1, time_step`, and closed once, after it. Inside, each former
!      `parallel do` becomes an `!$omp do` worksharing construct (which
!      reuses the already-live team), and the one truly serial bit per step
!      (updating the scalar E_val/E_vals/t_vals and the scalar observables
!      dipole_vals/population_vals + progress print) is wrapped in
!      `!$omp single`, whose implicit barrier keeps every thread correctly
!      synchronized on the shared state before/after it.
!
!   7) ANGULAR RECONSTRUCTION AND THE INTERACTION PHASE ARE FUSED.
!      Previously: psi_1(j,i) = sum_l Y_T(l,j)*A_mat(l,i) was written out in
!      full, then a second pass computed
!      psi_2(j,i) = exp(-i*E*cos(theta)*r*dt) * psi_1(j,i).
!      psi_1 required its own (L+1)x(N-1) complex array and its own sweep
!      over memory. Now both operations happen inside the same (j,i) loop
!      nest: the l_idx reduction is accumulated directly into psi_2(j,i),
!      and the interaction-phase factor is applied to that same accumulator
!      immediately afterwards, with no intermediate write/read of a psi_1
!      array. psi_1 is therefore eliminated completely.
!
!   8) tmp_glm IS ELIMINATED. The projection step used to accumulate into a
!      separate tmp_glm(l_idx,i) buffer and then divide into glm_tilde(i,l)
!      as a second statement. The reduction now accumulates directly into
!      glm_tilde(i,l_idx) and is normalized in place -- one array instead
!      of two, and one fewer array traversal per step.
!
!   9) THE INTERACTION-PHASE COEFFICIENT IS PRECOMPUTED. Every step the old
!      code recomputed cos_theta(j)*r(i)*dt (time-independent!) alongside
!      the genuinely time-dependent factor -E_val. This product is now
!      precomputed ONCE, before the time loop, into
!      phase_coeff(j,i) = cos_theta(j) * r(i) * dt
!      so that inside the ti-loop each step only has to do
!      phi = -E_val * phase_coeff(j,i)
!      i.e. one multiply instead of three, per (j,i) per step.
!
!  10) A NEW weighted_P(j,l) = weights(j) * P_table(j,l) TABLE IS
!      PRECOMPUTED, so the projection step's innermost loop no longer has
!      to multiply by weights(j) on every (j,l_idx,i) triple -- it was
!      constant across i and is now folded into the table once.
!
!  11) UNNECESSARY REPEATED CMPLX() CONVERSIONS ARE REMOVED FROM THE HOT
!      LOOPS. Fortran already allows direct real*complex arithmetic (the
!      real operand is promoted automatically), so wrapping every
!      Y_T(l_idx,j)*A_mat(l_idx,i) or weighted_P(j,l_idx)*psi_2(j,i) term in
!      cmplx(...) was pure overhead repeated inside the innermost reduction
!      loops. Likewise, the phase factor exp(cmplx(0.0d0, -phi, kind=8)) is
!      replaced with the equivalent cmplx(cos(phi), sin(phi), kind=8),
!      avoiding the more expensive complex EXP intrinsic on every grid
!      point.
!
!   Recommended compilation (see updated Makefile): compile with
!     -O3 -march=native -funroll-loops -fopenmp
!   and benchmark against the previous version. Control thread count at
!   run time with:
!     OMP_NUM_THREADS=<N> ./vector_time_evolution.out
!   Don't set N above l_max (or, for step 2, above L+1) -- extra threads
!   beyond the number of independent outer-loop iterations just add
!   scheduling overhead for no benefit.
! --------------------------------------------------------------------------

program main
    use parameters
    use timer_mod
    use legendre_stuff, only: a_legendre
    implicit none

    integer :: i, j, l_idx, ti, rec
    integer :: S_recl_size
    integer(kind=8), external :: factorial
    real(kind=8), external :: N_fact, C_fact
    character(len=256) :: colloc_file, state_file, smat_file, gl_file
    character(len=256) :: evo_data_file, evo_data_path
    character(len=8)   :: state_symb
    character(len=20) :: step_str, min_str, sec_str, tot_str

    ! [NOTE] E0_au, w0, cpp, dt, time_step, print_serial_prog, p_step, pi_au,
    ! evolving_atom, ... all now come from `use parameters` (parameters.f90) --
    ! nothing laser/system-related is hard-coded locally in this file anymore.
    ! Change the system in ONE place (parameters.f90: lambda_nm, I0, evolving_atom,
    ! confined, ...) and this program picks it up automatically on the next build.

    real(kind=8), dimension(L+1) :: roots, weights, cos_theta
    real(kind=8), dimension(N-1) :: x_glob, r, absorber, A_r
    real(kind=8), dimension(N-1, total_states) :: state_block

    ! --- Precomputed tables (built once, before the time loop) ---
    real(kind=8), dimension(l_max, L+1) :: Y_T          ! Y_T(l_idx, j)   -- l leading
    real(kind=8), dimension(L+1, l_max) :: P_table      ! P_table(j, l_idx) -- j leading

    ! weighted_P(j,l) = weights(j) * P_table(j,l)
    !
    ! Precomputed so the projection step's innermost loop does not have to
    ! multiply by weights(j) on every (j,l_idx,i) triple.
    real(kind=8), dimension(L+1, l_max) :: weighted_P

    real(kind=8), dimension(l_max) :: norm_factor        ! norm_factor(l_idx)

    ! phase_coeff(j,i) = cos_theta(j) * r(i) * dt
    !
    ! Time-independent part of the interaction phase. Inside the time loop
    ! only the genuinely time-dependent factor -E_val remains to be applied:
    !   phi = -E_val * phase_coeff(j,i)
    real(kind=8), dimension(L+1, N-1) :: phase_coeff

    ! --- matmul-facing arrays: l-index is the TRAILING dimension so slices
    !     taken at fixed l_idx are contiguous, cache-friendly matmul args ---
    complex(kind=8), dimension(N-1, N-1) :: S_l_matrix_c
    complex(kind=8), dimension(N-1, N-1, l_max+1) :: S_matrix_all
    complex(kind=8), dimension(N-1, l_max+1) :: init_glm, glm_tilde

    ! --- manual-reduction-loop-facing arrays: kept with l_idx leading,
    !     matching the innermost loop variable in those loops ---
    complex(kind=8), dimension(l_max+1, N-1) :: A_mat

    ! psi_1 has been eliminated: angular reconstruction and the interaction
    ! phase are now fused directly into psi_2 (see step 2 in the time loop).
    complex(kind=8), dimension(L+1, N-1) :: psi_2

    real(kind=8) :: exec_time, t_mid, E_val, phi
    real(kind=8), allocatable :: dipole_vals(:), population_vals(:)
    real(kind=8), allocatable :: t_vals(:), E_vals(:)

    integer :: pct, next_pct
    integer :: out_unit, io_stat

    ! time_step (how many steps to actually run) and time_step_max (the grid's
    ! full length) both come from parameters.f90 now -- see [Time evolution
    ! controls] there. Change `time_step` in parameters.f90 to run a shorter
    ! validation pass instead of the full pulse.
    allocate(dipole_vals(time_step), population_vals(time_step))
    allocate(t_vals(time_step), E_vals(time_step))

    ! ---------------------------------------------------------------
    ! 1) Radial grid: N-1 Gauss-Lobatto collocation points
    ! ---------------------------------------------------------------
    colloc_file = 'collocation_points/generator/Algo-3_Gauss_Lobatto_collocation_points.dat'
    if (.not. file_exists(colloc_file)) then
        colloc_file = '../collocation_points/generator/Algo-3_Gauss_Lobatto_collocation_points.dat'
    end if
    if (.not. file_exists(colloc_file)) then
        print *, 'ERROR: collocation file not found. Checked:'
        print *, '  - ', trim(colloc_file)
        stop 1
    end if

    open(unit=21, file=trim(colloc_file), status='old', action='read')
        read(21, *) x_glob
    close(21)

    ! ---------------------------------------------------------------
    ! 2) Angular grid: Gauss-Legendre nodes/weights for the theta grid
    ! ---------------------------------------------------------------
    gl_file = 'collocation_points/generator/Gauss_Legendre_collocation_points_and_weights.dat'
    if (.not. file_exists(gl_file)) then
        gl_file = '../collocation_points/generator/Gauss_Legendre_collocation_points_and_weights.dat'
    end if
    if (.not. file_exists(gl_file)) then
        print *, 'ERROR: Gauss-Legendre collocation/weight file not found. Checked:'
        print *, '  - ', trim(gl_file)
        stop 1
    end if

    open(unit=20, file=trim(gl_file), status='old', action='read')
        read(20, *) roots, weights
    close(20)

    r = Lmap * (1.0d0 + x_glob) / (1.0d0 - x_glob + alpha_map)
    cos_theta = roots

    state_file = 'GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin'
    if (.not. file_exists(state_file)) then
        state_file = '../GPSM_states_S-matrix/data_GPSM_states_S-matrix/GPSM-DSYEV_states.bin'
    end if
    if (.not. file_exists(state_file)) then
        print *, 'ERROR: GPSM state file not found. Checked:'
        print *, '  - ', trim(state_file)
        stop 1
    end if

    open(unit=22, file=trim(state_file), form='unformatted', access='stream', status='old')
        read(22) r
        read(22) state_block
    close(22)

    A_r = state_block(:, 1)
    init_glm = (0.0d0, 0.0d0)
    init_glm(:, 1) = cmplx(A_r, 0.0d0, kind=8)      ! l=0 partial wave -> first trailing index

    inquire(iolength=S_recl_size) S_l_matrix_c
    smat_file = 'GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin'
    if (.not. file_exists(smat_file)) then
        smat_file = '../GPSM_states_S-matrix/data_GPSM_states_S-matrix/S_matrices-DSYEV.bin'
    end if
    if (.not. file_exists(smat_file)) then
        print *, 'ERROR: S-matrix file not found. Checked:'
        print *, '  - ', trim(smat_file)
        stop 1
    end if

    open(unit=23, file=trim(smat_file), form='unformatted', access='direct', recl=S_recl_size, status='old')
    do rec = 1, l_max + 1
        read(23, rec=rec) S_l_matrix_c
        S_matrix_all(:, :, rec) = S_l_matrix_c            ! trailing l-index -> contiguous block per l
    end do
    close(23)

    do i = 1, N-1
        if (0.0d0 < r(i) .and. r(i) <= r0) then
            absorber(i) = 1.0d0
        else if (r(i) > r0 .and. r(i) < r_max) then
            absorber(i) = cos(pi_au * (r(i) - r0) / (2.0d0 * (r_max - r0)))**0.25d0
        else
            absorber(i) = 0.0d0
        end if
    end do

    ! ---------------------------------------------------------------
    ! Precompute EVERYTHING that does not change with ti, ONCE:
    !   - Y_T(l_idx, j)        : normalized Legendre values for the
    !                            angular reconstruction step (psi_2)
    !   - P_table(j, l_idx)    : raw Legendre values, same underlying
    !                            a_legendre() values as Y_T, just without
    !                            the N_fact factor, so we call a_legendre
    !                            only ONCE per (j,l_idx) pair total, not
    !                            once per (j,l_idx) PER TIME STEP.
    !   - weighted_P(j, l_idx) : weights(j) * P_table(j,l_idx), folding the
    !                            quadrature weight into the table so the
    !                            projection step's inner loop skips that
    !                            multiply.
    !   - norm_factor(l_idx)   : normalization used in the projection step,
    !                            independent of i, computed once per l.
    !   - phase_coeff(j, i)    : cos_theta(j)*r(i)*dt, the time-independent
    !                            part of the interaction phase.
    ! ---------------------------------------------------------------
    do l_idx = 1, l_max
        do j = 1, L+1
            P_table(j, l_idx) = a_legendre(l_idx-1, 0, roots(j))
            Y_T(l_idx, j) = N_fact(l_idx-1, 0) * P_table(j, l_idx)
            weighted_P(j, l_idx) = weights(j) * P_table(j, l_idx)
        end do
        norm_factor(l_idx) = N_fact(l_idx-1, 0) * C_fact(l_idx-1, 0)
    end do

    do j = 1, L+1
        do i = 1, N-1
            phase_coeff(j, i) = cos_theta(j) * r(i) * dt
        end do
    end do

    ! Human-readable state label (n+l, l) -> e.g. "1s"; fetched from parameters.f90::state_symbol
    state_symb = state_symbol(n_qn + l_qn, l_qn)

    print '(A)', '~~~~~~~~~~~: Time Evolution :~~~~~~~~~~'
    print '(A,A)', 'Evolving atom                   : ', trim(evolving_atom)
    print '(A,I0,A,I0,A,I0,A,A)', 'Evolving initial state (n,l,m)  : (', n_qn+l_qn, ', ', l_qn, ', ', m_qn, ') ~ ', trim(state_symb)
    print '(A,F10.2,A)', 'Wavelength (lambda_nm)          : ', lambda_nm, ' nm'
    print '(A,ES10.3,A)', 'Intensity (I0)                  : ', I0, ' W/cm^2'
    print '(A,I0)', 'Total time steps                : ', time_step
    print '(A)', 'Optimized OpenMP evolution enabled.'
    print *

    next_pct = 0
    call tick()

    ! ================================================================
    ! ================================================================
    !
    !     ONE PERSISTENT OPENMP TEAM AROUND THE ENTIRE TIME LOOP
    !
    ! The thread team is forked exactly once here (not four times per
    ! step). Every `!$omp do` below is a worksharing construct that reuses
    ! this already-live team; every `!$omp single` executes on one thread
    ! only, with its implicit barrier keeping the shared scalars/arrays
    ! correctly synchronized before the next worksharing construct reads
    ! them.
    !
    ! ================================================================
    ! ================================================================

    !$omp parallel default(shared) private(i, j, l_idx, ti, t_mid, phi)

    do ti = 1, time_step

        ! ------------------------------------------------------------
        ! Time quantities
        !
        ! Only one thread computes the scalar quantities. The implicit
        ! barrier at the end of SINGLE guarantees every thread sees the
        ! updated E_val before proceeding to step 1 below.
        !
        ! t_vals(ti) is the time AT the start of this step (matches
        ! Python's t[:time_step]); t_mid is used only for the interaction
        ! propagator, same as before. dt, E0_au, w0, cpp all come from
        ! parameters.f90.
        ! ------------------------------------------------------------

        !$omp single

        t_vals(ti) = dble(ti - 1) * dt
        t_mid = t_vals(ti) + 0.5d0 * dt

        E_val = E0_au * sin(w0 * t_mid) * sin(w0 * t_mid / (2.0d0 * cpp))**2

        ! Field value saved to file is E(t) at the step's start time, exactly
        ! like Python's saved column `E_field(t[:time_step])` (not the
        ! mid-point field).
        E_vals(ti) = E0_au * sin(w0 * t_vals(ti)) * sin(w0 * t_vals(ti) / (2.0d0 * cpp))**2

        !$omp end single

        ! ============================================================
        ! STEP 1
        !
        ! A_mat(l,:) = S(l) @ init_glm(:,l)
        !
        ! Both matmul args are contiguous (trailing l-index), so no hidden
        ! temp copies. Each l_idx reads its own S_matrix_all/init_glm slice
        ! and writes its own A_mat row -- fully independent across l_idx.
        ! ============================================================

        !$omp do schedule(static)

        do l_idx = 1, l_max
            A_mat(l_idx, :) = matmul(S_matrix_all(:, :, l_idx), init_glm(:, l_idx))
        end do

        !$omp end do

        ! ============================================================
        ! STEP 2 + INTERACTION PHASE (FUSED)
        !
        ! Original two passes:
        !   psi_1(j,i) = sum_l Y_T(l,j) * A_mat(l,i)
        !   psi_2(j,i) = exp(-i*E_val*cos(theta(j))*r(i)*dt) * psi_1(j,i)
        !
        ! are now fused into a single loop nest: the l_idx reduction
        ! accumulates directly into psi_2(j,i), and the interaction-phase
        ! factor is then applied in place. psi_1 is eliminated completely.
        !
        ! l_idx is the innermost (reduction) index; both Y_T and A_mat are
        ! laid out with l as the leading dimension for unit-stride access
        ! here. Parallelized over the outer indices (j,i): each (j,i) writes
        ! only its own psi_2 entry, so iterations are independent.
        !
        ! The phase factor uses cos/sin explicitly (instead of the complex
        ! EXP intrinsic) to avoid constructing a complex argument on every
        ! grid point, and no cmplx() wrapper is needed for the real*complex
        ! product Y_T(l_idx,j) * A_mat(l_idx,i) -- Fortran promotes the real
        ! operand automatically.
        ! ============================================================

        !$omp do collapse(2) schedule(static)

        do j = 1, L+1
            do i = 1, N-1

                psi_2(j, i) = (0.0d0, 0.0d0)

                do l_idx = 1, l_max
                    psi_2(j, i) = psi_2(j, i) + Y_T(l_idx, j) * A_mat(l_idx, i)
                end do

                phi = -E_val * phase_coeff(j, i)

                psi_2(j, i) = cmplx(cos(phi), sin(phi), kind=8) * psi_2(j, i)

            end do
        end do

        !$omp end do

        ! ============================================================
        ! STEP 3
        !
        ! Projection (tmp_glm eliminated -- accumulates directly into
        ! glm_tilde and normalizes in place):
        !
        ! glm_tilde(i,l) =
        !   [ sum_j weighted_P(j,l) * psi_2(j,i) ] / norm_factor(l)
        !
        ! weighted_P already has weights(j) folded in, so this inner loop
        ! does one multiply-add per j instead of two. j is the innermost
        ! (reduction) index; weighted_P and psi_2 are both laid out with j
        ! leading for unit-stride access here. No calls to a_legendre or
        ! N_fact/C_fact happen inside this loop anymore. Parallelized over
        ! the outer index l_idx: each l_idx writes only glm_tilde(:,l_idx),
        ! so iterations are independent.
        ! ============================================================

        !$omp do schedule(static)

        do l_idx = 1, l_max
            do i = 1, N-1

                glm_tilde(i, l_idx) = (0.0d0, 0.0d0)

                do j = 1, L+1
                    glm_tilde(i, l_idx) = glm_tilde(i, l_idx) + weighted_P(j, l_idx) * psi_2(j, i)
                end do

                glm_tilde(i, l_idx) = glm_tilde(i, l_idx) / norm_factor(l_idx)

            end do
        end do

        !$omp end do

        ! ============================================================
        ! STEP 4
        !
        ! init_glm(:,l) = S(l) @ glm_tilde(:,l) * absorber
        !
        ! Again both matmul args are contiguous trailing-l slices, parallel
        ! over l.
        ! ============================================================

        !$omp do schedule(static)

        do l_idx = 1, l_max
            init_glm(:, l_idx) = matmul(S_matrix_all(:, :, l_idx), glm_tilde(:, l_idx)) * absorber
        end do

        !$omp end do

        ! ============================================================
        ! SCALAR OBSERVABLES
        !
        ! Every worksharing construct above has an implicit barrier, so
        ! init_glm is complete before this point. Only one thread needs to
        ! compute these scalars and print progress.
        ! ============================================================

        !$omp single

        dipole_vals(ti) = dipole_scalar(r, init_glm)
        population_vals(ti) = population_scalar(init_glm)

        ! --- progress printing, mirrors parameters.py::print_serial_prog / p_step ---
        if (print_serial_prog) then
            pct = int(100.0d0 * dble(ti) / dble(time_step))
            if (pct >= next_pct) then
                print '(A,I8,A,A,I4,A,A)', 'Evolution step ', ti, ' : ', green, next_pct, '%', reset
                next_pct = next_pct + p_step
            end if
        end if

        !$omp end single

    end do

    !$omp end parallel

    call tock(exec_time)


    write(step_str, '(F12.5)') exec_time / dble(time_step)
    write(min_str, '(I0)') int(exec_time / 60.0d0)
    write(sec_str, '(F6.3)') mod(exec_time, 60.0d0)
    write(tot_str, '(I0)') time_step

    print '(A)', ''
    print '(A)', 'Average wall-time per step (eta_t) : ' // white // trim(adjustl(step_str)) // ' sec' // reset
    print '(A)', 'Total wall-time for all steps      : ' // green // trim(adjustl(min_str)) // ' min ' // &
        trim(adjustl(sec_str)) // ' sec' // reset
    print '(A)', 'Time steps used : ' // green // trim(adjustl(tot_str)) // reset
    ! ---------------------------------------------------------------
    ! Save t, E(t), d(t), Ps(t) -- same 4-column layout & header as
    ! vector_time_evolution.py (np.savetxt, fmt='%.16e').
    ! ---------------------------------------------------------------
    evo_data_file = 'VEvo_nopt='//trim(itoa(time_step))//'_'//trim(evolving_atom)//'('//trim(state_symb)//')_m='// &
                    trim(itoa(m_qn))//'_'//trim(SAE_model)//'_L='//trim(itoa(L))//'_kmax='//trim(itoa(kmax))// &
                    '_N='//trim(itoa(N))//'_rmax='//trim(itoa(int(r_max)))//'_Lmap='// &
                    trim(itoa(int(Lmap)))//'_dt='//trim(rtoa(dt))//'.dat'

    evo_data_path = 'Time_evolution_data/'//trim(evo_data_file)
    open(newunit=out_unit, file=trim(evo_data_path), status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
        ! Fall back to the current directory if 'Time_evolution_data/' doesn't exist here.
        evo_data_path = trim(evo_data_file)
        open(newunit=out_unit, file=trim(evo_data_path), status='replace', action='write', iostat=io_stat)
    end if

    if (io_stat /= 0) then
        print *, 'WARNING: could not open output file for writing: ', trim(evo_data_path)
    else
        write(out_unit, '(A)') 't(a.u.)       E(t)(a.u.)      d(t)(a.u.)      Ps(t)'
        do ti = 1, time_step
            write(out_unit, '(4(ES23.16E2,2X))') t_vals(ti), E_vals(ti), dipole_vals(ti), population_vals(ti)
        end do
        close(out_unit)
        print '(A,A)', "evo_data_file = ", orange // trim(evo_data_path) // reset
    end if


contains

    ! Small integer/real -> string helpers used for building the output filename.
    function itoa(i_val) result(s)
        integer, intent(in) :: i_val
        character(len=12) :: s
        write(s, '(I0)') i_val
    end function itoa

    function rtoa(r_val) result(s)
        real(kind=8), intent(in) :: r_val
        character(len=12) :: s
        write(s, '(F0.2)') r_val
    end function rtoa

    logical function file_exists(path)
        character(len=*), intent(in) :: path
        inquire(file=trim(path), exist=file_exists)
    end function file_exists

    ! NOTE: glm_arr is (N-1, l_max+1) -- radial index leading, l trailing --
    ! matching init_glm/glm_tilde used everywhere above.
    real(kind=8) function dipole_scalar(rp, glm_arr)
        real(kind=8), intent(in) :: rp(N-1)
        complex(kind=8), intent(in) :: glm_arr(N-1, l_max+1)
        integer :: l
        real(kind=8) :: sum_r, alpha_l
        sum_r = 0.0d0
        do l = 1, l_max
            alpha_l = sqrt(((dble(l))**2) / ((2.0d0 * dble(l) - 1.0d0) * (2.0d0 * dble(l) + 1.0d0)))
            sum_r = sum_r + alpha_l * sum(rp * real(conjg(glm_arr(:, l)) * glm_arr(:, l+1)))
        end do
        dipole_scalar = 2.0d0 * sum_r
    end function dipole_scalar

    real(kind=8) function population_scalar(glm_arr)
        complex(kind=8), intent(in) :: glm_arr(N-1, l_max+1)
        population_scalar = sum(abs(glm_arr)**2)
    end function population_scalar

end program main