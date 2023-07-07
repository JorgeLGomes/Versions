
# .SILENT:

include ../../configure/make.inc_gnu
#echo $(LIBDIR) is the lib directory

INC = ../../include

LINK.f = $(FC) $(FFLAGS) -I$(INC)

# Define main programs.

MAIN_TOPO = etatopo_3s.F
MAIN_TOPO2 = etatopo.F
MAIN_CORNERS = corners.f90

# Define other libraries needed to link.

#OTHER_LIB = $(LIBDIR)/libetautil.a 
OTHER_LIB = ../../../lib/libetautil.a 

# Define executables.

#EXE_TOPO = $(EXEDIR)/etatopo.exe
EXE_TOPO = ../../../exe/etatopo_3s.exe
EXE_TOPO2 = ../../../exe/etatopo.exe
EXE_CORNERS = ../../../exe/corners.exe

EXE = $(EXE_TOPO)
EXE2 = $(EXE_CORNERS)
EXE3 = $(EXE_TOPO2)

# Default target.

exe: $(EXE)
exe: $(EXE2)
exe: $(EXE3)

$(EXE_TOPO): $(MAIN_TOPO) $(OTHER_LIB)
	$(LINK.f) $(MAIN_TOPO) $(OTHER_LIB) -o $@
$(EXE_TOPO2): $(MAIN_TOPO2) $(OTHER_LIB)
	$(LINK.f) $(MAIN_TOPO2) $(OTHER_LIB) -o $@	
$(EXE_CORNERS): $(MAIN_CORNERS)
	$(FC) $(MAIN_CORNERS) -o $@
