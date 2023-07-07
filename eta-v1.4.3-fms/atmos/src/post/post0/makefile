################################################################################
include	../../configure/make.inc
EXE=	../../../exe/post0.x

LDFLAGS  =
#
# Non-threaded object files
#
OBJS=	HORIZ.o POST0.o PTETAE.o VRT_TSLP.o
#
# Includes
#
INCLUDES=  parmeta
#
# Common Blocks
#
COMMS= 	

DEPS= $(COMMS) $(INCLUDES)

.SUFFIXES:	.f90 .F .f .o

.f.o:
	$(FC) $(FFLAGS) -c $<	

.f90.o:
	$(FC) $(FFLAGS) -c $<	

$(EXE):	$(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

clean:	
	/bin/rm -f  *.lst *.o *.i
