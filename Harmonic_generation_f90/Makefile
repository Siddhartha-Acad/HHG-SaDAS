# Makefile for Harmonic_generation_f90
# gl_nodes  : Gauss–Legendre nodes/collocation points
# glo_nodes : Gauss–Lobatto nodes/collocation points

CYAN := \033[1;36m
YELLOW := \033[1;33m
GREEN := \033[1;32m
RED := \033[1;31m
RESET := \033[0m

# --------------------------------------------------------------------
# [PERFORMANCE] Previously plain `gfortran` (~-O0), a 5-20x slowdown
# vs -O3 for numerically heavy code. Applied to every compile/link step.
#
# -O3/-march=native/-funroll-loops : aggressive optimization for this CPU
# -fopenmp : drives the !$omp parallel/do directives in
#            vector_time_evolution.f90 (independent per-l/j/i slices)
#
# -ffast-math kept OFF (changes results via relaxed IEEE rounding) --
# see FAST_FFLAGS below for a separate benchmarking-only variant.
# --------------------------------------------------------------------
FFLAGS := -O3 -march=native -funroll-loops -fopenmp

# Optional aggressive floating-point flags, for benchmarking only.
# Always numerically validate against a plain -O3 run before trusting
# results built with this.
FAST_FFLAGS := $(FFLAGS) -ffast-math

# --------------------------------------------------------------------
# BLAS
#
# -fexternal-blas routes vector_time_evolution.f90's large MATMUL calls
# (the S(l) applies) through real BLAS instead of gfortran's naive loop.
# -fblas-matmul-limit=0 removes the small-matrix cutoff so it applies
# even at the (N-1)x(N-1) sizes used here.
#
# Scoped to only vector_time_evolution.out's link step, not global
# FFLAGS, so the OpenBLAS dependency isn't pulled into the smaller
# generator programs that don't need it.
#
# Requires: sudo apt install libopenblas-dev
# --------------------------------------------------------------------
VEC_BLAS_FFLAGS := -fexternal-blas -fblas-matmul-limit=0
VEC_BLAS_LIBS   := -lopenblas

# --------------------------------------------------------------------
# [OPENMP THREADS] Default OMP_NUM_THREADS=2 for `make run_...`
# (override: `make run_vector_time_evolution OMP_NUM_THREADS=8`, or
# set it yourself when running the binary directly).
#
# OPENBLAS_NUM_THREADS=1 / MKL_NUM_THREADS=1 pin BLAS to 1 thread so
# it doesn't spawn its own team inside the already-parallel OpenMP
# region (oversubscription). OMP_DYNAMIC=FALSE stops the runtime from
# silently shrinking the team size.
# --------------------------------------------------------------------
OMP_NUM_THREADS ?= 2
export OMP_NUM_THREADS

OPENBLAS_NUM_THREADS ?= 1
export OPENBLAS_NUM_THREADS

MKL_NUM_THREADS ?= 1
export MKL_NUM_THREADS

OMP_DYNAMIC ?= FALSE
export OMP_DYNAMIC

COLLOC_dPATH := collocation_points/generator
INPUT_dPATH  := GPSM_states_S-matrix/data_GPSM_states_S-matrix
INPUT_PATH   := GPSM_states_S-matrix
CHECK_PATH   := correctness_check

# --------------------------------------------------------------------
# Stamp file marking a successful vector_time_evolution.out run.
# vector_time_evolution.f90 writes many Time_evolution_data/*.dat
# files rather than one fixed name, so a glob can't be used directly
# as a Make target with reliable timestamps. Touching this stamp
# after a successful run gives Make one file whose mtime it can
# compare against vector_time_evolution.out, so HHG_analysis (and
# run_vector_time_evolution itself) only reruns the simulation when
# the binary is actually newer than the last completed run.
# --------------------------------------------------------------------
TIME_EVOL_STAMP := Time_evolution_data/.stamp

.PHONY : all clean clobber glo_nodes gl_nodes GPSM_states S_matrices vector_time_evolution run_vector_time_evolution HHG_analysis check_GPSM check_Split_operator

all : gl_nodes glo_nodes GPSM_states S_matrices run_vector_time_evolution HHG_analysis

gl_nodes    : $(COLLOC_dPATH)/Gauss_Legendre_collocation_points_and_weights.dat
glo_nodes   : $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_collocation_points.dat
GPSM_states : glo_nodes $(INPUT_dPATH)/GPSM-DSYEV_states.bin
S_matrices  : glo_nodes $(INPUT_dPATH)/S_matrices-DSYEV.bin


$(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_collocation_points.dat : $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto.f90 functions.o parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                    Algo-3_Gauss_Lobatto.f90                      |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)gfortran $(FFLAGS) -cpp $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto.f90 functions.o measure_time.o -o $(COLLOC_dPATH)/Gauss_Lobatto.out$(RESET)"
	@gfortran $(FFLAGS) -cpp -I. $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto.f90 functions.o measure_time.o -o $(COLLOC_dPATH)/Gauss_Lobatto.out
	@echo "$(YELLOW)cd $(COLLOC_dPATH) && ./Gauss_Lobatto.out$(RESET)"
	@cd $(COLLOC_dPATH) && ./Gauss_Lobatto.out

$(COLLOC_dPATH)/Gauss_Legendre_collocation_points_and_weights.dat : $(COLLOC_dPATH)/Gauss_Legendre.f90 parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|           Gauss_Legendre.f90 :: Golub-Welsch Algorithm           |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)gfortran $(FFLAGS) $(COLLOC_dPATH)/Gauss_Legendre.f90 -llapack -o $(COLLOC_dPATH)/Gauss_Legendre.out$(RESET)"
	@gfortran $(FFLAGS) $(COLLOC_dPATH)/Gauss_Legendre.f90 measure_time.o -llapack -o $(COLLOC_dPATH)/Gauss_Legendre.out
	@echo "$(YELLOW)cd $(COLLOC_dPATH) && ./Gauss_Legendre.out$(RESET)"
	@cd $(COLLOC_dPATH) && ./Gauss_Legendre.out

$(INPUT_dPATH)/GPSM-DSYEV_states.bin : $(INPUT_PATH)/GPSM-DSYEV.f90 parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                         GPSM-DSYEV.f90                           |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(INPUT_PATH) && gfortran $(FFLAGS) -I.. GPSM-DSYEV.f90 ../measure_time.o -llapack -o GPSM_gen.out$(RESET)"
	@cd $(INPUT_PATH) && gfortran $(FFLAGS) -I.. GPSM-DSYEV.f90 ../measure_time.o -llapack -o GPSM_gen.out
	@echo "$(YELLOW)cd $(INPUT_PATH) && ./GPSM_gen.out$(RESET)"
	@cd $(INPUT_PATH) && ./GPSM_gen.out

vector_time_evolution : vector_time_evolution.out

# --------------------------------------------------------------------
# vector_time_evolution.out is the only target linked against BLAS
# (VEC_BLAS_FFLAGS / VEC_BLAS_LIBS), since it's the only program with
# MATMUL calls large enough to benefit -- see the BLAS comment above.
#
# It also reads GPSM-DSYEV_states.bin and S_matrices-DSYEV.bin at
# runtime, so both are now listed as prerequisites. This means a bare
# `make run_vector_time_evolution` on a fresh clone will transitively
# build glo_nodes -> GPSM_states -> S_matrices -> vector_time_evolution
# in the right order, instead of segfaulting on missing/garbage input.
# --------------------------------------------------------------------
vector_time_evolution.out : vector_time_evolution.f90 parameters.o functions.o measure_time.o \
                            $(INPUT_dPATH)/GPSM-DSYEV_states.bin $(INPUT_dPATH)/S_matrices-DSYEV.bin
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                 vector_time_evolution.f90                         |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)gfortran $(FFLAGS) $(VEC_BLAS_FFLAGS) -I. vector_time_evolution.f90 parameters.o functions.o measure_time.o $(VEC_BLAS_LIBS) -o vector_time_evolution.out$(RESET)"
	@gfortran $(FFLAGS) $(VEC_BLAS_FFLAGS) -I. vector_time_evolution.f90 parameters.o functions.o measure_time.o $(VEC_BLAS_LIBS) -o vector_time_evolution.out


# --------------------------------------------------------------------
# [STACK SIZE] Large automatic arrays / MATMUL temporaries (sized off
# the (N-1)x(N-1) S-matrices) live on the stack, and the default 8MB
# limit is too small, causing a silent SIGSEGV. `ulimit -s unlimited`
# fixes this for the main thread; run via `bash -c` since ulimit is a
# shell builtin, not an external command.
#
# Workaround, not a permanent fix -- durable fix is making those
# arrays ALLOCATABLE (heap) in vector_time_evolution.f90.
#
# --------------------------------------------------------------------
# [UP-TO-DATE CHECK] The real target here is the stamp file, not the
# .out binary itself, so Make can tell whether the *last completed
# run* is newer than the current binary. If vector_time_evolution.out
# hasn't changed since the last run, this recipe is skipped entirely
# -- `make run_vector_time_evolution` becomes a no-op, and downstream
# `make HHG_analysis` goes straight to plotting.
# --------------------------------------------------------------------
run_vector_time_evolution : $(TIME_EVOL_STAMP)

$(TIME_EVOL_STAMP) : vector_time_evolution.out
	@mkdir -p Time_evolution_data
	@echo "$(YELLOW)OMP_NUM_THREADS=$(OMP_NUM_THREADS) ulimit -s unlimited ./vector_time_evolution.out$(RESET)"
	@bash -c "ulimit -s unlimited && ./vector_time_evolution.out"
	@touch $(TIME_EVOL_STAMP)


# --------------------------------------------------------------------
# Depends on the stamp file (real timestamp), not the phony
# run_vector_time_evolution alias, so `make HHG_analysis` only
# reruns the simulation when vector_time_evolution.out is newer than
# the last completed run -- otherwise it goes straight to plotting
# the existing Time_evolution_data/*.dat files.
# --------------------------------------------------------------------
HHG_analysis : HHG_analysis/HHG_spectra.py $(TIME_EVOL_STAMP)
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                        HHG_spectra.py                            |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)python3 HHG_analysis/HHG_spectra.py$(RESET)"
	@python3 HHG_analysis/HHG_spectra.py



$(INPUT_dPATH)/S_matrices-DSYEV.bin : $(INPUT_PATH)/S-matrix-DSYEV.f90 parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                       S-matrix-DSYEV.f90                         |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(INPUT_PATH) && gfortran $(FFLAGS) -cpp -I.. S-matrix-DSYEV.f90 ../measure_time.o -llapack -o S-matrix_gen.out$(RESET)"
	@cd $(INPUT_PATH) && gfortran $(FFLAGS) -cpp -I.. S-matrix-DSYEV.f90 ../measure_time.o -llapack -o S-matrix_gen.out
	@echo "$(YELLOW)cd $(INPUT_PATH) && ./S-matrix_gen.out$(RESET)"
	@cd $(INPUT_PATH) && ./S-matrix_gen.out

$(INPUT_dPATH)/Eigenstates-DSYEV.bin : $(INPUT_PATH)/S-matrix-DSYEV.f90 parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                 S-matrix-DSYEV.f90 -DSAVE_STATES                 |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(INPUT_PATH) && gfortran $(FFLAGS) -cpp -DSAVE_STATES -I.. S-matrix-DSYEV.f90 ../measure_time.o -llapack -o A_states_gen.out$(RESET)"
	@cd $(INPUT_PATH) && gfortran $(FFLAGS) -cpp -DSAVE_STATES -I.. S-matrix-DSYEV.f90 ../measure_time.o -llapack -o A_states_gen.out
	@echo "$(YELLOW)cd $(INPUT_PATH) && ./A_states_gen.out$(RESET)"
	@cd $(INPUT_PATH) && ./A_states_gen.out

check_GPSM  : GPSM_states
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                         check_GPSM.f90                           |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(CHECK_PATH) && gfortran $(FFLAGS) -I.. check_GPSM.f90 ../measure_time.o -o check_GPSM.out$(RESET)"
	@cd $(CHECK_PATH) && gfortran $(FFLAGS) -I.. check_GPSM.f90 ../measure_time.o -o check_GPSM.out
	@echo "$(YELLOW)cd $(CHECK_PATH) && ./check_GPSM.out$(RESET)"
	@cd $(CHECK_PATH) && ./check_GPSM.out

check_Split_operator : S_matrices $(INPUT_dPATH)/Eigenstates-DSYEV.bin
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                    check_Split_operator.f90                      |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(CHECK_PATH) && gfortran $(FFLAGS) -I.. check_Split_operator.f90 ../measure_time.o -o check_Split_operator.out$(RESET)"
	@cd $(CHECK_PATH) && gfortran $(FFLAGS) -I.. check_Split_operator.f90 ../measure_time.o -o check_Split_operator.out
	@echo "$(YELLOW)cd $(CHECK_PATH) && ./check_Split_operator.out$(RESET)"
	@cd $(CHECK_PATH) && ./check_Split_operator.out


parameters.o : parameters.f90
	@echo "$(YELLOW)gfortran $(FFLAGS) -c parameters.f90$(RESET)"
	@gfortran $(FFLAGS) -c parameters.f90

functions.o : functions.f90
	@echo "$(YELLOW)gfortran $(FFLAGS) -c functions.f90$(RESET)"
	@gfortran $(FFLAGS) -c functions.f90

measure_time.o : measure_time.f90
	@echo "$(YELLOW)gfortran $(FFLAGS) -c measure_time.f90$(RESET)"
	@gfortran $(FFLAGS) -c measure_time.f90


clean :
	@echo "$(RED)rm -f $(COLLOC_dPATH)/Gauss_Lobatto.out$(RESET)"
	@rm -f $(COLLOC_dPATH)/Gauss_Lobatto.out
	@echo "$(RED)rm -f $(COLLOC_dPATH)/Gauss_Legendre.out$(RESET)"
	@rm -f $(COLLOC_dPATH)/Gauss_Legendre.out
	@echo "$(RED)rm -f $(INPUT_PATH)/GPSM_gen.out$(RESET)"
	@rm -f $(INPUT_PATH)/GPSM_gen.out
	@echo "$(RED)rm -f $(INPUT_PATH)/S-matrix_gen.out$(RESET)"
	@rm -f $(INPUT_PATH)/S-matrix_gen.out
	@echo "$(RED)rm -f $(INPUT_PATH)/A_states_gen.out$(RESET)"
	@rm -f $(INPUT_PATH)/A_states_gen.out
	@echo "$(RED)rm -f $(CHECK_PATH)/check_GPSM.out$(RESET)"
	@rm -f $(CHECK_PATH)/check_GPSM.out
	@echo "$(RED)rm -f $(CHECK_PATH)/check_Split_operator.out$(RESET)"
	@rm -f $(CHECK_PATH)/check_Split_operator.out
	@echo "$(RED)rm -f vector_time_evolution.out$(RESET)"
	@rm -f vector_time_evolution.out

clobber : clean
	@echo "$(RED)rm -f *.o *.mod$(RESET)"
	@rm -f *.o *.mod
	@echo "$(RED)rm -f $(COLLOC_dPATH)/*.dat$(RESET)"
	@rm -f $(COLLOC_dPATH)/*.dat
	@echo "$(RED)rm -f $(INPUT_dPATH)/*.bin$(RESET)"
	@rm -f $(INPUT_dPATH)/*.bin
	@echo "$(RED)rm -f Time_evolution_data/*$(RESET)"
	@rm -f Time_evolution_data/*
