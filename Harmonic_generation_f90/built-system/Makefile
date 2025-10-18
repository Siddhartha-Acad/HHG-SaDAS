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
# %   : To create any .o file, first look for a .f90 file with the same base name.
# $<  : if target is main.o, $< = main.f90