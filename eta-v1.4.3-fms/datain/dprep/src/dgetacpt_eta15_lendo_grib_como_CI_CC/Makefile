
.SILENT:

RM = rm -f

# Define compiler flags from machine dependent make.inc.

include ../configure/make.inc

# Define main routine.

MAIN = dgeta2model.f

# Define executable.

EXE = $(EXEDIR)/dgcpteta15.exe

OTHERLIB = ../../../libraries/w3lib
# Default target

$(EXE): $(MAIN)
	$(FC) $(FFLAGS) $(MAIN) $(OTHERLIB) -o $@
	$(RM) *.o

exe: $(EXE)


# Define sst interp main routine.

MAIN2 = interpola_sst.f90

# Define sst interp executable.

EXE2 = $(EXEDIR)/interpola_sst.x

$(EXE2): $(MAIN2)
	$(FC) $(MAIN2) -o $(EXE2)
	$(RM) *.o

exe: $(EXE2)





