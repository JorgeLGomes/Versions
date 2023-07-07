#!/bin/sh

#
#     Generate a list of object files that corresponds to the
#     list of Fortran ( .f ) files in the current directory
#

# Include compiler options ------------------------------------
include ../../eta/src/configure/make.inc
# by Lamosa PAD/CPTEC/2007 ------------------------------------
list=`ls *.f | grep '.f' `
for routine in $list
do
name=`echo $routine | cut -f1 -d.`
	if ${FC} -c -o $name.o $name.f
	then
	ar r ../bufrlib $name.o
	rm -f $name.o
	fi
done
