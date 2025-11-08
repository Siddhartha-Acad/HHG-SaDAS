# Makefile for Harmonic_generation_f90

COLLOC_dPATH := ./collocation_points/generator
INPUT_dPATH  := ./GPSM_states_S-matrix/data_GPSM_states_S-matrix
INPUT_PATH   := ./GPSM_states_S-matrix
CHECK_PATH   := ./correctness_ckeck

.PHONY : all clear clobber colloc_pt GPSM_states S_matrices

all : colloc_pt GPSM_states S_matrices

colloc_pt   : $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_collocation_points.dat
GPSM_states : colloc_pt $(INPUT_dPATH)/GPSM-DSYEV_states.bin
S_matrices  : colloc_pt $(INPUT_dPATH)/S_matrices-DSYEV.bin
check_GPSM  : GPSM_states
	cd $(CHECK_PATH) && gfortran -I.. check_GPSM.f90 -o check_GPSM.out
	cd $(CHECK_PATH) && ./check_GPSM.out

$(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_collocation_points.dat : $(COLLOC_dPATH)/Algo-3.f90 functions.o
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "|                           Algo-3.f90                             |"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	gfortran -cpp -I. $(COLLOC_dPATH)/Algo-3.f90 functions.o -o $(COLLOC_dPATH)/Algo-3.out
	cd $(COLLOC_dPATH) && ./Algo-3.out

$(INPUT_dPATH)/GPSM-DSYEV_states.bin : $(INPUT_PATH)/GPSM-DSYEV.f90 parameters.o
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "|                         GPSM-DSYEV.f90                           |"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	cd $(INPUT_PATH) && gfortran -I.. GPSM-DSYEV.f90 -llapack -lblas -o GPSM_gen.out
	cd $(INPUT_PATH) && ./GPSM_gen.out

$(INPUT_dPATH)/S_matrices-DSYEV.bin : $(INPUT_PATH)/S-matrix-DSYEV.f90 parameters.o
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "|                       S-matrix-DSYEV.f90                         |"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	cd $(INPUT_PATH) && gfortran -cpp -I.. S-matrix-DSYEV.f90 -llapack -lblas -o S-matrix_gen.out
	cd $(INPUT_PATH) && ./S-matrix_gen.out

parameters.o : parameters.f90
	gfortran -c parameters.f90

functions.o : functions.f90
	gfortran -c functions.f90

clean :
	rm -f *.o *.mod
	rm -f $(COLLOC_dPATH)/Algo-3.out
	rm -f $(INPUT_PATH)/GPSM_gen.out
	rm -f $(INPUT_PATH)/S-matrix_gen.out
	rm -f $(CHECK_PATH)/check_GPSM.out

clobber : clean
	rm -f $(COLLOC_dPATH)/*.dat
	rm -f $(INPUT_dPATH)/*.bin

