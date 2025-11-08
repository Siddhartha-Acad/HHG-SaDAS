# Makefile for Harmonic_generation_f90

COLLOC_dPATH := ./collocation_points/generator
INPUT_dPATH  := ./GPSM_states_S-matrix/data_GPSM_states_S-matrix
GPSM_PATH    := ./GPSM_states_S-matrix


.PHONY : clear clobber colloc_pt GPSM_states

$(COLLOC_dPATH)/Algo-3.out : $(COLLOC_dPATH)/Algo-3.f90 functions.o
	gfortran -I. $(COLLOC_dPATH)/Algo-3.f90 functions.o -o $(COLLOC_dPATH)/Algo-3.out

parameters.o : parameters.f90
	gfortran -c parameters.f90

functions.o : functions.f90
	gfortran -c functions.f90

colloc_pt : $(COLLOC_dPATH)/Algo-3.out
	cd $(COLLOC_dPATH) && ./Algo-3.out


$(GPSM_PATH)/GPSM_gen.out : $(GPSM_PATH)/GPSM-DSYEV.f90 parameters.o
	cd $(GPSM_PATH) && gfortran -I.. GPSM-DSYEV.f90 ../parameters.o -llapack -lblas -o GPSM_gen.out

$(INPUT_dPATH)/GPSM-DSYEV_states.bin : $(GPSM_PATH)/GPSM_gen.out
	cd $(GPSM_PATH) && ./GPSM_gen.out

GPSM_states : $(INPUT_dPATH)/GPSM-DSYEV_states.bin


clean :
	rm -f *.o *.mod
	rm -f $(COLLOC_dPATH)/Algo-3.out
	rm -f $(GPSM_PATH)/GPSM_gen.out

clobber : clean
	rm -f $(COLLOC_dPATH)/*.dat
	rm -f $(INPUT_dPATH)/*.bin

