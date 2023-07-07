
.SILENT:

RM = rm -f

# Define compiler flags from machine dependent make.inc.
include ../../configure/make.inc
COMPILE.f = $(FC) -c 

# Define main routine.
CMD=	../../../exe/reform_2d3d.x	
SRC=	variables.f90 constants.f90 TQadj.f90   thetae.f90   reform_2d3d.f90 
OBJ=	variables.o constants.o TQadj.o   thetae.o   reform_2d3d.o 
################################################################################
.SUFFIXES: .f90 .o
.f90.o:
	echo "Compiling $(<:.f90=.f90)"
	$(COMPILE.f) $(<:.f90=.f90) 
$(CMD):	$(OBJ)
	$(FC) -o $@ $(OBJ)
clean:
	/bin/rm -f  *.lst *.o *.mod
exe:	$(CMD)
