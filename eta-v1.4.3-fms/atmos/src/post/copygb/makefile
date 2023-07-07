SHELL=	/bin/sh
include ../../configure/make.inc
CMD=	../../../exe/copygb.x
LOPTS=	
INCS=	
OBJS=	copygb.o  
LIBA= ../../../../libraries/w3lib 
LIBC= ../../../../libraries/iplib
LIBS= $(LIBA) $(LIBC) 
################################################################################
$(CMD):	$(OBJS)
	$(FC) -o $(CMD) $(OBJS) $(LIBS)
$(OBJS):	$(INCS)
copygb.o:	copygb.f
	$(FC) -c $(FFLAGS) $*.f
clean:	
	-rm -f $(OBJS) *.lst *.i
