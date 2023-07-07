# Start of script compip.sh          M.Farley  96-10-22
#
echo Compile and place in library iplib subroutine $1.f
#
# Example:  Compile subroutine q9yi32.f and put it into library iplib
#
#           compip.sh q9yi32
#
rm    $1.o               # delete subroutine old object
if f90 -c -O nofastint $1.o $1.f # compile subroutine
then
ar r ../iplib $1.o
rm    $1.o               # delete subroutine new object
rm    $1.l               # delete listing file
#
echo Subroutine $1.f was compiled and placed in library iplib
else
echo Compile error, subroutine $1.f was not placed in library iplib
fi
#
# End of script compip.sh
