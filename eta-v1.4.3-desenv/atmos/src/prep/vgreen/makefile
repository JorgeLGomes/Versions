################################################################################
# 
#     Makefile for VGREEN - 
#
#	Updated for version 1.3.0 of the Eta model Eta
#
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#     The following macros will be of interest:
#
#         TARGET   - name of the executable
#         FC       - name of Fortran compiler
#         CPP      - name of CPP
#         CPPFLAGS - CPP flags
#         OPTS     - compiler code optimizations
#         LIST     - source listing
#         SMP      - threading
#         TRAPS    - runtime traps for floating point exceptions
#         PROFILE  - source code profiling ( -pg )
#         DEBUG    - -g
#         INCLUDE  - search location for include files
#         W3LIB    - w3lib
#         MPI      - dummy MPI library
#         SEARCH   - library search location
#
#################################################################################

.SUFFIXES: .f90 .o
# INC = ../include_mpi
INC=./
include ../../configure/make.inc

#
# Define the name of the executable
#
TARGET = ../../../exe/vgreen.x

COMPILE.f = $(FC) -c $(FFLAGS) -I$(INC)
PREPROC.F = $(CPP) -P $(CFLAGS)

CPP      = /lib/cpp -P
CPPFLAGS = 
OPTS= -O1
SMP      = 
TRAPS    =
PROFILE  =
DEBUG=
INCLUDE = ./
W3LIB= 
# MPI=../../libraries/dummyMPI/libmpi.a
# SEARCH   = -L../../libraries/dummyMPI
#
# Assemble Options
#
# FFLAGS should be defined in make.inc
# FFLAGS   = $(OPTS) $(LIST) $(TRAPS) $(PROFILE) $(DEBUG) $(INCLUDE)
FFLAGST  = $(FFLAGS) $(SMP)
LDFLAGS  = $(SMP) $(PROFILE) $(FFLAGS)
#LIBS     = $(MPI) $(W3LIB) -lmpi
LIBS     = 
#
# Non-threaded object files
#
OBJS=	MODULE_F77KINDS.o MODULE_PARMCONF.o MODULE_CONSTANTS.o MODULE_DIAGNOSTIC.o\
	GEOGRAPHIC.o PUTVEG.o VFRAC.o   \

#
# Threaded object files
#
OBJST=  

#
# Includes
#
INCLUDES=  
#
# Common Blocks
#
COMMS=  

DEPS= $(COMMS) $(INCLUDES)

.SUFFIXES: .F .o

.f90.o:
#	echo "USING CPP STUFF"
#	$(PREPROC.F) $< > $*.f
#	mv $(<:.F=.i)  $(<:.F=.f)
	echo "Compiling $(<:.f90=.f90)"
	$(COMPILE.f) $(<:.f90=.f90)
#	rm $(<:.F=.f)

$(TARGET):	$(OBJS) $(OBJST)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(OBJST) $(LIBS)
#	$(FC) $(LDFLAGS) $(LIBS) -o $@ $(OBJS) $(OBJST)

# $(OBJS):	$(DEPS)
#	$(COMPILE.f) $<

# $(OBJST):	$(DEPS)
#	$(FC) $(FFLAGST) -c $<

clean:	
	/bin/rm -f  *.i *.lst *.o *.mod

exe: $(TARGET)
