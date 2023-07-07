
.SILENT:

RM = rm -f

# Define compiler flags from machine dependent make.inc.

include ../configure/make.inc

# Define main routine.

MAIN = dgeta2model.f

# Define executable.

EXE = $(EXEDIR)/dgcpteta40.exe

OTHERLIB = ../../../../libraries/w3lib

# Default target

$(EXE): $(MAIN)
	$(FC) $(FFLAGS) $(MAIN) $(OTHERLIB) -o $@
	$(RM) *.o

exe: $(EXE)

