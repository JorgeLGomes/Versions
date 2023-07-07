################################################################################
# 
#     Makefile for Eta Coordinate Quilt
#
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#################################################################################

F77=f90

#
include ../../configure/make.inc
#
# Define the name of the executable
#
TARGET = ../../../exe/quilt.x
#
# CPP, Compiler, and Linker Options
#
# CPP      = /lib/cpp -P
MPI      = ../../../../libraries/dummyMPI/libmpi.a
#
# Assemble Options
#
LDFLAGS  = 
LIBS     = $(MPI)
#
#
# Non-threaded object files
#
OBJS= QUILT.o SLP.o SLPSIG.o SLPSIGSPLINE.o PARA_RANGE.o MPI_FIRST.o \
	 MPI_LAST.o COLLECT.o DIST.o \
	 UPDATE.o
#
# Includes
#
INCLUDES=  parmeta parmsoil parmbuf
#
# Common Blocks
#
# COMMS=	PARA.comm

# DEPS= $(COMMS) $(INCLUDES)

.SUFFIXES:	.F .f .o

$(TARGET):	$(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

# $(OBJS):	$(DEPS)
# $(OBJS):
#	$(FC) $(FFLAGS) -c $<

clean:	
	/bin/rm -f  $(TARGET) *.lst *.o lm map *.i
