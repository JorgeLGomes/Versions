
# .SILENT:

include ../../configure/make.inc
# echo $(LIBDIR) is the lib directory

INC = ../../include

LINK.f = $(FC) $(FFLAGS) -I$(INC)

# Define main programs.

MAIN_VEG = vegmsk.f

# Define other libraries needed to link.

OTHER_LIB = 

# Define executables.

EXE_VEG = $(EXEDIR)/vegmsk.exe

EXE = $(EXE_VEG)

# Default target.

exe: $(EXE)

$(EXE_VEG): $(MAIN_VEG) $(OTHER_LIB)
	$(LINK.f) $(MAIN_VEG) $(OTHER_LIB) -o $@

clean:
	/bin/rm -f *i. *.lst *.o *.i
