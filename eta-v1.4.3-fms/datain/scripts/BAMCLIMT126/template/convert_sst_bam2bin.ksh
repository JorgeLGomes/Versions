#!/bin/ksh -x
if (($#<1)) ; then
   echo 'Parameter: yyyymmddhh'
   exit
fi
Run_Date=$1
export Dir_scr=DIRROOT
export Dir_home=`dirname ${Dir_scr}`
export Dir_datain=`dirname ${Dir_home}`
export ModelDrive=`basename ${Dir_scr}`
export Dir_wrk=${Dir_datain}/sst/${ModelDrive}/${Run_Date}
export Dir_ETAwrk_SST=${Dir_datain}/sst/ETAwrk/${ModelDrive}/${Run_Date}
export Dir_wrkt=${Dir_ETAwrk_SST}/WRK
export Dir_util=${Dir_datain}/util

if [ ! -d ${Dir_ETAwrk_SST} ] ; then
  mkdir -p  ${Dir_ETAwrk_SST}
fi
if [ ! -d ${Dir_wrkt} ] ; then
  mkdir -p  ${Dir_wrkt}
fi

typeset -Z4 TNTimes
cdo -f nc4 import_binary ${Dir_wrk}/INPUTSST${Run_Date:0:8}.G00192.ctl  ${Dir_wrkt}/INPUTSST${Run_Date:0:8}.G00192.nc
cdo select,name=sstp -addc,273.15 ${Dir_wrkt}/INPUTSST${Run_Date:0:8}.G00192.nc  ${Dir_wrkt}/INPUTSST${Run_Date}.G00192_TempK.nc
cdo -sellonlatbox,0,360,90,-90 ${Dir_wrkt}/INPUTSST${Run_Date}.G00192_TempK.nc ${Dir_wrkt}/INPUTSST${Run_Date}.G00192_TempK_0-360_90_-90.nc
cdo -f grb -copy ${Dir_wrkt}/INPUTSST${Run_Date}.G00192_TempK_0-360_90_-90.nc ${Dir_wrkt}/INPUTSST${Run_Date}.G00192_TempK_0-360_90_-90.grib
TNTimes=`cdo -ntime ${Dir_wrkt}/INPUTSST${Run_Date}.G00192_TempK_0-360_90_-90.nc`
echo ${TNTimes}
echo ${TNTimes} >${Dir_ETAwrk_SST}/NTimes.txt

i=1
while (($i<=TNTimes));do
  ${Dir_util}/wgrib -d $i ${Dir_wrkt}/INPUTSST${Run_Date}.G00192_TempK_0-360_90_-90.grib -nh -bin -append -o ${Dir_wrkt}/wfile_0-360_-90_90.bin
  let i=$i+1
done

mv ${Dir_wrkt}/wfile_0-360_-90_90.bin ${Dir_ETAwrk_SST}/${ModelDrive}_${Run_Date}.bin 
rm -Rf ${Dir_wrkt}
  exit
