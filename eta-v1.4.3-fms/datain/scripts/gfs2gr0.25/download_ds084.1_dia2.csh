#!/bin/csh
#################################################################
# Csh Script to retrieve 372 online Data files of 'ds084.1',
# total 125.42G. This script uses 'wget' to download data.
#
# Highlight this script by Select All, Copy and Paste it into a file;
# make the file executable and run it on command line.
#
# You need pass in your password as a parameter to execute
# this script; or you can set an environment variable RDAPSWD
# if your Operating System supports it.
#
# Contact rpconroy@ucar.edu (Riley Conroy) for further assistance.
#################################################################

set rundate=$2
set pswd = $1
if(x$pswd == x && `env | grep RDAPSWD` != '') then
 set pswd = $RDAPSWD
endif
if(x$pswd == x) then
 echo
 echo Usage: $0 YourPassword
 echo
 exit 1
endif
set v = `wget -V |grep 'GNU Wget ' | cut -d ' ' -f 3`
set a = `echo $v | cut -d '.' -f 1`
set b = `echo $v | cut -d '.' -f 2`
if(100 * $a + $b > 109) then
 set opt = 'wget --no-check-certificate'
else
 set opt = 'wget'
endif
set opt1 = '-O Authentication.log --save-cookies auth.rda_ucar_edu --post-data'
set opt2 = "email=diego.chagas@inpe.br&passwd=$pswd&action=login"
$opt $opt1="$opt2" https://rda.ucar.edu/cgi-bin/login
set opt1 = "-N --load-cookies auth.rda_ucar_edu"
set opt2 = "$opt $opt1 https://rda.ucar.edu/data/ds084.1/"
set filelist = ( \
  2021/$rundate/gfs.0p25.${rundate}00.f000.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f006.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f012.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f018.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f024.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f030.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f036.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f042.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f048.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f054.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f060.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f066.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f072.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f078.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f084.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f090.grib2 \
  2021/$rundate/gfs.0p25.${rundate}00.f096.grib2 \
 )
while($#filelist > 0)
 set syscmd = "$opt2$filelist[1]"
 echo "$syscmd ..."
 $syscmd
 shift filelist
end

