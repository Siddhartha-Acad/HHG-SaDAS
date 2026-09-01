# Makefile for Harmonic_generation_f90
# gl_nodes  : Gauss–Legendre nodes/collocation points
# glo_nodes : Gauss–Lobatto nodes/collocation points

CYAN := \033[1;36m
YELLOW := \033[1;33m
GREEN := \033[1;32m
RED := \033[1;31m
RESET := \033[0m

# --------------------------------------------------------------------
# [PERFORMANCE] The previous version of this Makefile compiled every
# target with plain `gfortran` and no optimization flags at all, which
# is equivalent to -O0. That alone can be a 5-20x slowdown vs -O3 for
# numerically heavy code like this. FFLAGS below is now applied to every
# compile/link step.
#
# -O3              : aggressive optimization
# -march=native    : use the instruction set of the current CPU
# -funroll-loops   : unroll suitable hot loops
# -fopenmp         : enable OpenMP (drives the "!$omp parallel"/"!$omp do"
#                    directives around the per-l S-matrix apply loops and
#                    the fused reconstruction/projection loops in
#                    vector_time_evolution.f90 -- each l/j/i slice is
#                    independent, so this scales with core count)
#
# Keep -ffast-math OFF here because this is a physics calculation and
# -ffast-math can change results (it relaxes IEEE rounding/NaN/Inf
# handling). A separate FAST_FFLAGS below lets you benchmark that
# variant explicitly, on its own binary, without touching production
# builds.
# --------------------------------------------------------------------
FFLAGS := -O3 -march=native -funroll-loops -fopenmp

# Optional aggressive floating-point flags, for benchmarking only.
# Always numerically validate against a plain -O3 run before trusting
# results built with this.
FAST_FFLAGS := $(FFLAGS) -ffast-math

# --------------------------------------------------------------------
# BLAS
#
# vector_time_evolution.f90 contains large MATMUL operations (the S(l)
# applies in steps 1 and 4). With gfortran, -fexternal-blas routes
# suitable MATMUL expressions through a real BLAS implementation instead
# of gfortran's own naive internal loop, which can be a large speed-up.
# -fblas-matmul-limit=0 removes the small-matrix cutoff so this kicks in
# even for the (N-1)x(N-1) sizes used here.
#
# These flags/libs are scoped to ONLY the vector_time_evolution.out link
# step below (VEC_BLAS_FFLAGS / VEC_BLAS_LIBS), rather than added to the
# global FFLAGS -- this avoids pulling the OpenBLAS dependency into the
# smaller generator programs (Gauss_Lobatto, Gauss_Legendre, GPSM_gen,
# S-matrix_gen, ...) that don't need it.
#
# Requires:
#     sudo apt install libopenblas-dev
# --------------------------------------------------------------------
VEC_BLAS_FFLAGS := -fexternal-blas -fblas-matmul-limit=0
VEC_BLAS_LIBS   := -lopenblas

# --------------------------------------------------------------------
# [OPENMP THREADS] Default thread count used when running
# vector_time_evolution.out through `make run_vector_time_evolution`.
# `?=` means: use 2 unless OMP_NUM_THREADS is already set in the
# environment or passed on the command line, e.g.:
#   make run_vector_time_evolution OMP_NUM_THREADS=8
# `export` makes it visible to the ./vector_time_evolution.out child
# process. Running the binary directly (not via `make run_...`) is
# unaffected by this default -- set OMP_NUM_THREADS yourself in that
# case, e.g. OMP_NUM_THREADS=8 ./vector_time_evolution.out
#
# [AVOIDING THREAD OVERSUBSCRIPTION] vector_time_evolution.f90 already
# parallelizes over l (and over j,i) with OpenMP. If MATMUL is also
# routed through a threaded BLAS (see VEC_BLAS_LIBS above), that BLAS
# would by default try to spawn its OWN thread team INSIDE the already-
# parallel OpenMP region, oversubscribing the machine (N_omp x N_blas
# threads competing for N_cpu cores) and making things slower, not
# faster. OPENBLAS_NUM_THREADS=1 / MKL_NUM_THREADS=1 pin the BLAS calls
# made from inside the OpenMP region to a single thread each, so all
# parallelism comes from the outer OpenMP loops as intended.
# OMP_DYNAMIC=FALSE stops the OpenMP runtime from silently shrinking the
# team size, so OMP_NUM_THREADS is respected exactly for benchmarking.
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

.PHONY : all clear clobber glo_nodes gl_nodes GPSM_states S_matrices vector_time_evolution run_vector_time_evolution HHG_analysis \

all : gl_nodes glo_nodes GPSM_states S_matrices vector_time_evolution

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
# --------------------------------------------------------------------
vector_time_evolution.out : vector_time_evolution.f90 parameters.o functions.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                 vector_time_evolution.f90                         |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)gfortran $(FFLAGS) $(VEC_BLAS_FFLAGS) -I. vector_time_evolution.f90 parameters.o functions.o measure_time.o $(VEC_BLAS_LIBS) -o vector_time_evolution.out$(RESET)"
	@gfortran $(FFLAGS) $(VEC_BLAS_FFLAGS) -I. vector_time_evolution.f90 parameters.o functions.o measure_time.o $(VEC_BLAS_LIBS) -o vector_time_evolution.out

run_vector_time_evolution : vector_time_evolution
	@echo "$(YELLOW)OMP_NUM_THREADS=$(OMP_NUM_THREADS) ./vector_time_evolution.out$(RESET)"
	@./vector_time_evolution.out

HHG_analysis : HHG_analysis/HHG_spectra.py
	@echo "$(YELLOW)HHG analysis script is ready in HHG_analysis/$(RESET)"

HHG_analysis/HHG_spectra.py :
	@mkdir -p HHG_analysis
	@touch HHG_analysis/HHG_spectra.py

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
