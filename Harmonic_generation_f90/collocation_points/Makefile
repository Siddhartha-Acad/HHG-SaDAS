# Algo-3.f90: to generate GL collocation points.

ifeq ($(OS), Windows_NT)
    RM = del
    EXE = Algo-3.exe
else
    RM = rm
    EXE = ./Algo-3.exe
endif


.PHONY : colloc_pt clear


colloc_pt : Algo-3.exe
	cd ./generator && $(EXE) && $(RM) Algo-3.exe

Algo-3.exe : ./generator/Algo-3.f90 newton_raphson.o functions.o
	gfortran -I.. ./generator/Algo-3.f90 ./newton_raphson.o ./functions.o -o ./generator/Algo-3.exe

functions.o : ../functions.f90
	gfortran -J.. -c ../functions.f90

newton_raphson.o : newton_raphson.f90 functions.o
	gfortran -I.. -c newton_raphson.f90

clear :
ifeq ($(OS), Windows_NT)
	del ..\legendre_stuff.mod
	del .\generator\*.dat
	del *.o
else
	rm ../legendre_stuff.mod
	rm ./generator/*.dat
	rm *.o
endif