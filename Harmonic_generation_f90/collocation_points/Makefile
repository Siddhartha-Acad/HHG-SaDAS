# Algo-3.f90: to generate Gauss-Lobatto collocation points.

.PHONY : colloc_pt clear clear_dat

./generator/Algo-3.exe : ./generator/Algo-3.f90 functions.o parameters.o
	gfortran -I.. ./generator/Algo-3.f90 functions.o parameters.o -o ./generator/Algo-3.exe

parameters.o : ../parameters.f90
	gfortran -J.. -c ../parameters.f90

functions.o : ../functions.f90
	gfortran -J.. -c ../functions.f90


colloc_pt : ./generator/Algo-3.exe
	cd ./generator && ./Algo-3.exe

clear :
	rm -f *.o
	rm -f ./generator/Algo-3.exe
	rm -f ../legendre_stuff.mod ../parameters.mod

clear_dat :
	rm -f ./generator/*.dat
