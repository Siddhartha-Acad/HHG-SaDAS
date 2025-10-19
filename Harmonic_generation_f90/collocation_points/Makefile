# Algo-3.f90: to generate GL collocation points.

.PHONY : colloc_pt clear clear_all


colloc_pt : Algo-3.exe
	./generator/Algo-3.exe

Algo-3.exe : ./generator/Algo-3.f90 functions.o newton_raphson.o
	gfortran -I.. ./generator/Algo-3.f90 ./functions.o ./newton_raphson.o -o ./generator/Algo-3.exe

functions.o : ../functions.f90
	gfortran -J.. -c ../functions.f90

newton_raphson.o : newton_raphson.f90
	gfortran -I.. -c newton_raphson.f90

clear :
	rm ./generator/Algo-3.exe

clear_all :
	rm ../legendre_stuff.mod
	rm ./generator/Algo-3.exe
	rm *.dat
	rm *.o
