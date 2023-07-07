# Makefile 
#
# by Lamosa PAD/CPTEC/2007

# Include compiler options ------------------------------------
include ../configure/make.inc
list=`ls *.f | grep '.f'`

all: clean
	@${FC} -c *.f
	@ar r ../bufrlib *.o
	@rm -f *.o
	@echo "ei_bufrlib library complete..."

clean:
	@rm -f ../bufrlib

