################################################################################
# 
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#     The following macros will be of interest:
#
#         TARGET   - name of the executable
#         FC       - name of Fortran compiler
#         OPTS     - compiler code optimizations
#         LIST     - source listing
#         SMP      - threading
#         TRAPS    - runtime traps for floating point exceptions
#         PROFILE  - source code profiling ( -pg )
#         DEBUG    - -g
#         SEARCH   - library search location
#
################################################################################
#

.SUFFIXES: .f90 .o
#
# Define the name of the executable
#
include ../configure/make.inc
TARGET = $(EXEDIR)/dgCFS.exe

FFLAGS = 
CPPFLAGS = 
OPTS= -O1
SMP      = 
TRAPS    =
PROFILE  =
DEBUG=
INCLUDE = ./
#
# Assemble Options
#
# FFLAGS should be defined in make.inc
# FFLAGS   = $(OPTS) $(LIST) $(TRAPS) $(PROFILE) $(DEBUG) $(INCLUDE)
FFLAGST  = $(FFLAGS) $(SMP)
LDFLAGS  = $(SMP) $(PROFILE) $(FFLAGS)
OTHERLIB = 
#
# Non-threaded object files
#
OBJS=	module_estab.o module_constants.o conv2model.o \
	es_ini.o rh_fix.o TQadj.o Input2ModelData.o
#
# Common Blocks
#
.SUFFIXES: .f90 .o 

.f90.o:
	echo "Compiling $(<:.f90=.f90)"
	$(COMPILE.f) $(<:.f90=.f90)

$(TARGET):	$(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)   $(OTHERLIB)
clean:	
	/bin/rm -f  *.lst *.o *.mod

exe: $(TARGET)
