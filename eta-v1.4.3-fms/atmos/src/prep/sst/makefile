################################################################################
# 
#     Makefile for Eta Coordinate Model - Dummy MPI Version
#                                          Tuccillo Dec 1999
# 	
#	Updated for internal quilting version of the model
#					   Pyle July 2000
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
TARGET = ../../../exe/sst.x

COMPILE.f = $(FC) -c -O0 $(FFLAGS) -I$(INC)
PREPROC.F = $(CPP) -P $(CFLAGS)

CPP      = /lib/cpp -P
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
FFLAGST  = $(FFLAGS) $(SMP)
LDFLAGS  = $(SMP) $(PROFILE) $(FFLAGS)
LIBS     = 
#
# Non-threaded object files
#
OBJS=	MODULE_F77KINDS.o MODULE_PARMCONF.o MODULE_CONSTANTS.o MODULE_DIAGNOSTIC.o \
	SST_FILL.o Geographic.o READ_SST.o CED_IJ.o SSTG2EGRD.o  \

#
# Includes
#
INCLUDES=  
#

DEPS= $(COMMS) $(INCLUDES)

.SUFFIXES: .f90 .o

.f90.o:
#	echo "USING CPP STUFF"
	echo "Compiling $(<:.f90=.f90)"
	$(COMPILE.f) $(<:.f90=.f90)

$(TARGET):	$(OBJS) $(OBJST)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(OBJST) $(LIBS)
#	$(FC) $(LDFLAGS) $(LIBS) -o $@ $(OBJS) $(OBJST)

clean:	
	/bin/rm -f  *.i *.lst *.o *.mod

exe: $(TARGET)
