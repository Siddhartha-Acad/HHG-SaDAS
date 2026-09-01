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
# If you have OpenBLAS or MKL installed, you can get a further large
# speed-up on the matmul() calls in vector_time_evolution.f90 by letting
# gfortran route its intrinsic MATMUL through a real BLAS implementation
# instead of its own naive internal loop. To enable that, uncomment the
# two BLAS_* lines below (requires e.g. `sudo apt install libopenblas-dev`).
#
# -fopenmp enables the "!$omp parallel do" directives around the per-l
# S-matrix apply loops in vector_time_evolution.f90 (each l is independent,
# so this scales with core count). Control thread count at run time with:
#   OMP_NUM_THREADS=8 ./vector_time_evolution.out
# --------------------------------------------------------------------
FFLAGS := -O3 -march=native -funroll-loops -fopenmp

# BLAS_FFLAGS := -fexternal-blas -fblas-matmul-limit=0
# BLAS_LIBS   := -lopenblas
BLAS_FFLAGS :=
BLAS_LIBS   :=

COLLOC_dPATH := collocation_points/generator
INPUT_dPATH  := GPSM_states_S-matrix/data_GPSM_states_S-matrix
INPUT_PATH   := GPSM_states_S-matrix
CHECK_PATH   := correctness_ckeck

.PHONY : all clear clobber glo_nodes gl_nodes GPSM_states S_matrices vector_time_evolution run_vector_time_evolution HHG_analysis

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
	@echo "$(CYAN)|           Gauss_Legendre.f90 :: Golub–Welsch Algorithm           |$(RESET)"
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

vector_time_evolution.out : vector_time_evolution.f90 parameters.o functions.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                 vector_time_evolution.f90                         |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)gfortran $(FFLAGS) $(BLAS_FFLAGS) -I. vector_time_evolution.f90 parameters.o functions.o measure_time.o $(BLAS_LIBS) -o vector_time_evolution.out$(RESET)"
	@gfortran $(FFLAGS) $(BLAS_FFLAGS) -I. vector_time_evolution.f90 parameters.o functions.o measure_time.o $(BLAS_LIBS) -o vector_time_evolution.out

run_vector_time_evolution : vector_time_evolution
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

clobber : clean
	@echo "$(RED)rm -f *.o *.mod$(RESET)"
	@rm -f *.o *.mod
	@echo "$(RED)rm -f $(COLLOC_dPATH)/*.dat$(RESET)"
	@rm -f $(COLLOC_dPATH)/*.dat
	@echo "$(RED)rm -f $(INPUT_dPATH)/*.bin$(RESET)"
	@rm -f $(INPUT_dPATH)/*.bin