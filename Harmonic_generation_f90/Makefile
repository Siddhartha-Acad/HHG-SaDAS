# Makefile for Harmonic_generation_f90

colloc_path := ./collocation_points/generator

.PHONY : colloc_pt clear clobber

$(colloc_path)/Algo-3.out : $(colloc_path)/Algo-3.f90 functions.o
	gfortran -I. $(colloc_path)/Algo-3.f90 functions.o -o $(colloc_path)/Algo-3.out

parameters.o : parameters.f90
	gfortran -c parameters.f90

functions.o : functions.f90
	gfortran -c functions.f90


colloc_pt : $(colloc_path)/Algo-3.out
	cd $(colloc_path) && ./Algo-3.out

clear :
	rm -f *.o
	rm -f $(colloc_path)/Algo-3.out
	rm -f legendre_stuff.mod
	rm -f parameters.mod parameters.o

clobber :
	rm -f $(colloc_path)/*.dat

