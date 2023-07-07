#
#   File     :  makeprofile
#   Written  :  30 September 2010, Jorge Luis Gomes
#   Revisions:  none
#
#   This makefile creates executables of the eta profile post processor.
#
#
#
#   Set source and object files.
#
#
RM = rm -f 
include ../../configure/make.inc
CMD=    ../../../exe/profile.x
SRC=    parmeta.f90 profile.f90
#

$(CMD):	$(SRC) 
	$(FC) $(SRC) -o $(@)
	$(RM) *.o *.mod

exe:	$(CMD)

#   Option to remove object files (*.o) after make.
#
clean: 
#   rm -f *.[omjl]

