SHELL=	/bin/sh
include	../../configure/make.inc
CMD=	../../../exe/staids.x
LOPTS=	
INCS=	
OBJS=	STALST.o
####################################################
$(CMD): $(OBJS)
	$(FC) -o $(CMD) $(OBJS)
$(OBJS):	$(INCS)
STALST.o:	STALST.f
	$(FC) -c $(FFLAGS) $*.f
clean:	
	-rm -f $(OBJS) *.lst *.i
