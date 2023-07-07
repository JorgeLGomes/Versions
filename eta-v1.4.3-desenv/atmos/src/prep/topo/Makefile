
# .SILENT:

include ../../configure/make.inc_initbc
# echo $(LIBDIR) is the lib directory

INC = ../../include

LINK.f = $(FC) $(FFLAGS) -I$(INC)

# Define main programs.

MAIN_TOPO = etatopo.F

# Define other libraries needed to link.

OTHER_LIB = $(LIBDIR)/libetautil.a 

# Define executables.

EXE_TOPO = $(EXEDIR)/etatopo.exe

EXE = $(EXE_TOPO)

# Default target.

exe: $(EXE)

$(EXE_TOPO): $(MAIN_TOPO) $(OTHER_LIB)
	$(LINK.f) $(MAIN_TOPO) $(OTHER_LIB) -o $@
