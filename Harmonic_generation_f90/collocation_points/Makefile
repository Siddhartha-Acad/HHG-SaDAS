# Algo-3.f90: to generate Gauss-Lobatto collocation points.

.PHONY : colloc_pt clear clear_dat

./generator/Algo-3.exe : ./generator/Algo-3.f90 newton_raphson.o functions.o
	gfortran -I.. ./generator/Algo-3.f90 ./newton_raphson.o ./functions.o -o ./generator/Algo-3.exe

functions.o ../legendre_stuff.mod: ../functions.f90
	gfortran -J.. -c ../functions.f90

newton_raphson.o : newton_raphson.f90 ../legendre_stuff.mod
	gfortran -I.. -c newton_raphson.f90

colloc_pt : ./generator/Algo-3.exe
	cd ./generator && ./Algo-3.exe

clear :
	rm -f *.o
	rm -f ./generator/Algo-3.exe
	rm -f ../legendre_stuff.mod

clear_dat :
	rm -f ./generator/*.dat
