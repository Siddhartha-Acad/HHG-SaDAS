#                Makefile
# +----------------------------------------+
# | target: dependencies                   |
# | <TAB>  command(s) to make target       |
# +----------------------------------------+

FC = gfortran
OBJECTS = main.o sub1.o sub2.o

output.txt : main.exe
	.\main.exe > output.txt

main.exe : $(OBJECTS)
	$(FC) $(OBJECTS) -o main.exe

%.o : %.f90
	$(FC) -c $<

# FC  : Fortran Compiler
# %   : for any .o file it requires a .f90 file.
# $<  : whatever appear in dependency side.