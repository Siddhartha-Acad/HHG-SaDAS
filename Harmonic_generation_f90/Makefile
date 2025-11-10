# Makefile for Harmonic_generation_f90

CYAN := \033[1;36m
YELLOW := \033[1;33m
GREEN := \033[1;32m
RED := \033[1;31m
RESET := \033[0m

COLLOC_dPATH := ./collocation_points/generator
INPUT_dPATH  := ./GPSM_states_S-matrix/data_GPSM_states_S-matrix
INPUT_PATH   := ./GPSM_states_S-matrix
CHECK_PATH   := ./correctness_ckeck

.PHONY : all clear clobber colloc_pt GPSM_states S_matrices

all : colloc_pt GPSM_states S_matrices

colloc_pt   : $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_collocation_points.dat
GPSM_states : colloc_pt $(INPUT_dPATH)/GPSM-DSYEV_states.bin
S_matrices  : colloc_pt $(INPUT_dPATH)/S_matrices-DSYEV.bin


$(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_collocation_points.dat : $(COLLOC_dPATH)/Algo-3.f90 functions.o parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                           Algo-3.f90                             |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)gfortran -cpp -I../.. $(COLLOC_dPATH)/Algo-3.f90 functions.o measure_time.o -o $(COLLOC_dPATH)/Algo-3.out$(RESET)"
	@gfortran -cpp -I../.. $(COLLOC_dPATH)/Algo-3.f90 functions.o measure_time.o -o $(COLLOC_dPATH)/Algo-3.out
	@echo "$(YELLOW)cd $(COLLOC_dPATH) && ./Algo-3.out$(RESET)"
	@cd $(COLLOC_dPATH) && ./Algo-3.out

$(INPUT_dPATH)/GPSM-DSYEV_states.bin : $(INPUT_PATH)/GPSM-DSYEV.f90 parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                         GPSM-DSYEV.f90                           |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(INPUT_PATH) && gfortran -I.. GPSM-DSYEV.f90 ../measure_time.o -llapack -lblas -o GPSM_gen.out$(RESET)"
	@cd $(INPUT_PATH) && gfortran -I.. GPSM-DSYEV.f90 ../measure_time.o -llapack -lblas -o GPSM_gen.out
	@echo "$(YELLOW)cd $(INPUT_PATH) && ./GPSM_gen.out$(RESET)"
	@cd $(INPUT_PATH) && ./GPSM_gen.out

$(INPUT_dPATH)/S_matrices-DSYEV.bin : $(INPUT_PATH)/S-matrix-DSYEV.f90 parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                       S-matrix-DSYEV.f90                         |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(INPUT_PATH) && gfortran -cpp -I.. S-matrix-DSYEV.f90 ../measure_time.o -llapack -lblas -o S-matrix_gen.out$(RESET)"
	@cd $(INPUT_PATH) && gfortran -cpp -I.. S-matrix-DSYEV.f90 ../measure_time.o -llapack -lblas -o S-matrix_gen.out
	@echo "$(YELLOW)cd $(INPUT_PATH) && ./S-matrix_gen.out$(RESET)"
	@cd $(INPUT_PATH) && ./S-matrix_gen.out

$(INPUT_dPATH)/Eigenstates-DSYEV.bin : $(INPUT_PATH)/S-matrix-DSYEV.f90 parameters.o measure_time.o
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                 S-matrix-DSYEV.f90 -DSAVE_STATES                 |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(INPUT_PATH) && gfortran -cpp -DSAVE_STATES -I.. S-matrix-DSYEV.f90 ../measure_time.o -llapack -lblas -o A_states_gen.out$(RESET)"
	@cd $(INPUT_PATH) && gfortran -cpp -DSAVE_STATES -I.. S-matrix-DSYEV.f90 ../measure_time.o -llapack -lblas -o A_states_gen.out
	@echo "$(YELLOW)cd $(INPUT_PATH) && ./A_states_gen.out$(RESET)"
	@cd $(INPUT_PATH) && ./A_states_gen.out

check_GPSM  : GPSM_states
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                         check_GPSM.f90                           |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(CHECK_PATH) && gfortran -I.. check_GPSM.f90 ../measure_time.o -o check_GPSM.out$(RESET)"
	@cd $(CHECK_PATH) && gfortran -I.. check_GPSM.f90 ../measure_time.o -o check_GPSM.out
	@echo "$(YELLOW)cd $(CHECK_PATH) && ./check_GPSM.out$(RESET)"
	@cd $(CHECK_PATH) && ./check_GPSM.out

check_Split_operator : S_matrices $(INPUT_dPATH)/Eigenstates-DSYEV.bin
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                    check_Split_operator.f90                      |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)cd $(CHECK_PATH) && gfortran -I.. check_Split_operator.f90 ../measure_time.o -o check_Split_operator.out$(RESET)"
	@cd $(CHECK_PATH) && gfortran -I.. check_Split_operator.f90 ../measure_time.o -o check_Split_operator.out
	@echo "$(YELLOW)cd $(CHECK_PATH) && ./check_Split_operator.out$(RESET)"
	@cd $(CHECK_PATH) && ./check_Split_operator.out

parameters.o : parameters.f90
	@echo "$(YELLOW)gfortran -c parameters.f90$(RESET)"
	@gfortran -c parameters.f90

functions.o : functions.f90
	@echo "$(YELLOW)gfortran -c functions.f90$(RESET)"
	@gfortran -c functions.f90

measure_time.o : measure_time.f90
	@echo "$(YELLOW)gfortran -c measure_time.f90$(RESET)"
	@gfortran -c measure_time.f90


clean :
	@echo "$(RED)rm -f $(COLLOC_dPATH)/Algo-3.out$(RESET)"
	@rm -f $(COLLOC_dPATH)/Algo-3.out
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

