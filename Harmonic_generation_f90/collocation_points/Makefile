# Algo-3.f90: to generate GL collocation points.

.PHONY : colloc_pt clear clear_dat

colloc_pt : ./generator/Algo-3.exe
	cd ./generator && ./Algo-3.exe

./generator/Algo-3.exe : ./generator/Algo-3.f90 newton_raphson.o functions.o
	gfortran -I.. ./generator/Algo-3.f90 ./newton_raphson.o ./functions.o -o ./generator/Algo-3.exe

functions.o ../legendre_stuff.mod: ../functions.f90
	gfortran -J.. -c ../functions.f90

newton_raphson.o : newton_raphson.f90 ../legendre_stuff.mod
	gfortran -I.. -c newton_raphson.f90

clear :
	rm *.o
	rm ./generator/Algo-3.exe
	rm ../legendre_stuff.mod

clear_dat :
	rm ./generator/*.dat
