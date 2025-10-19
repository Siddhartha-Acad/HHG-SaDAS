# Algo-3.f90: to generate GL collocation points.

RM = del

.PHONY : colloc_pt clear


colloc_pt : Algo-3.exe
	./generator/Algo-3.exe
	$(RM) .\generator\Algo-3.exe

Algo-3.exe : ./generator/Algo-3.f90 newton_raphson.o functions.o
	gfortran -I.. ./generator/Algo-3.f90 ./newton_raphson.o ./functions.o -o ./generator/Algo-3.exe

functions.o : ../functions.f90
	gfortran -J.. -c ../functions.f90

newton_raphson.o : newton_raphson.f90 functions.o
	gfortran -I.. -c newton_raphson.f90

clear :
	$(RM) ..\legendre_stuff.mod
	$(RM) .\generator\Algo-3.exe
	$(RM) *.dat
	$(RM) *.o
