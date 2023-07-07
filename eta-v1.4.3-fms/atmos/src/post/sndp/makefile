#
SHELL=/bin/sh
include ../../configure/make.inc
F77=$(FC)
#
SRCS=	BFRHDR.f BFRIZE.f CALHEL.f CALWXT.f SMIX.f SNDPST.f VAP.f \
        WETBLB.f CALVIS.f 

OBJS=	BFRHDR.o BFRIZE.o CALHEL.o CALWXT.o SMIX.o SNDPST.o VAP.o \
        WETBLB.o CALVIS.o 

# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling
#
LDFLAGS = 
CMD =		../../../exe/sndp.x


LIBS =  ../../../../libraries/w3lib ../../../../libraries/bufrlib 

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

# Make the profiled version of the command and call it a.out.prof
#
$(CMD).prof:	$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	-rm -f $(OBJS) *.lst *.i

clobber:	clean
	-rm -f $(CMD) $(CMD).prof

void:	clobber
	-rm -f $(SRCS) makefile
