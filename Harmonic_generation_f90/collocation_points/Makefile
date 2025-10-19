# Algo-3.f90: to generate GL collocation points.

run : Algo-3.exe
	./Algo-3.exe

Algo-3.exe : Algo-3.f90 functions.o newton_raphson.o
	gfortran -I.. Algo-3.f90 functions.o newton_raphson.o -o Algo-3.exe

functions.o : ../functions.f90
	gfortran -J.. -c ../functions.f90

newton_raphson.o : newton_raphson.f90
	gfortran -I.. -c newton_raphson.f90
