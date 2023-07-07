# Makefile iplib 

# by Lamosa PAD/CPTEC/2007

# Include copmpilers options ------------------------------------
include ../configure/make.inc

all: clean
	@${FC}  -c *.f
	@ar r ../iplib *.o
	@rm -f *.o
	@echo "iplib library complete..."

clean:
	@rm -f ../iplib
