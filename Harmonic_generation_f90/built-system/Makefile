#                Makefile
# +----------------------------------------+
# | target: dependencies                   |
# | <TAB>  command(s) to make target       |
# +----------------------------------------+

OBJECTS = main.o sub1.o sub2.o

output.txt : main.exe
	.\main.exe > output.txt

main.exe : $(OBJECTS)
	gfortran $(OBJECTS) -o main.exe

%.o : %.f90
	gfortran -c $<

# %   means: 'for any .o file it requires a .f90 file'.
# $<  means whatever appear in dependency side.