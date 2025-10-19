# Algo-3.f90: to generate GL collocation points.

run : Algo-3.exe
	./Algo-3.exe

Algo-3.exe : Algo-3.f90 newton_raphson.o functions.o
	gfortran -I.. Algo-3.f90 newton_raphson.o functions.o -o Algo-3.exe

functions.o : ../functions.f90
	gfortran -J.. -c ../functions.f90

newton_raphson.o : newton_raphson.f90 functions.o
	gfortran -I.. -c newton_raphson.f90
