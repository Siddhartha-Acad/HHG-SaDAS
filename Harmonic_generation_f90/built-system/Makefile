#                Makefile
# +----------------------------------------+
# | target: dependencies                   |
# | <TAB>  command(s) to make target       |
# +----------------------------------------+


output.txt : main.exe
	.\main.exe > output.txt

main.exe : main.o sub1.o sub2.o
	gfortran main.o sub1.o sub2.o -o main.exe

%.o : %.f90
	gfortran -c $<

# %   means: 'for any .o file it requires a .f90 file'.
# $<  means whatever appear in dependency side.