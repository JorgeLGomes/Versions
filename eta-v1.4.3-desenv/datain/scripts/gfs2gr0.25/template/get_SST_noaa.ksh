#!/bin/ksh -x
#
#
# Verifica se foi passado algum parametro
hostname
if (($#<1)) ; then
export Run_Date=`date "+%Y%m%d"`00
else
export Run_Date=${1}00
fi

# VARIAVEIS
export Dir_scr=DIRROOT
export Dir_home=`dirname ${Dir_scr}`
export Dir_datain=`dirname ${Dir_home}`
export Dir_util=${Dir_datain}/util
export IntSST=24

cd ${Dir_scr}

datadir=pub/data/nccf/com/nsst/prod
server=https://nomads.ncep.noaa.gov
dir_grb2=nsst

datesst=`${Dir_util}/caldate.3.0 ${Run_Date} - 1dy 'yyyymmdd'`
datesst2=`${Dir_util}/caldate.3.0 ${Run_Date} - 1dy 'yyyymmddhh'`

export Dir_SST=${Dir_datain}/sst/noaa/${datesst}
export Dir_ETAwrk_SST=${Dir_datain}/sst/ETAwrk/noaa/${datesst2}
if [ ! -d ${Dir_SST} ] ; then
  mkdir -p  ${Dir_SST}
fi
if [ ! -d ${Dir_ETAwrk_SST} ] ; then
  mkdir -p  ${Dir_ETAwrk_SST}
fi


cd ${Dir_SST}
ifile=rtgssthr_grb_0.083
dir_grb2=nsst.${datesst} 
if [[ ! -s ${Dir_SST}/${ifile}_${datesst}.grib2 ]] ; then
   echo $server/${datadir}/${dir_grb2}/${ifile}.grib2 
   wget --no-check-certificate -q -c $server/${datadir}/${dir_grb2}/${ifile}.grib2 -O ${ifile}_${datesst}.grib2 
fi
if  [ -s ${Dir_SST}/${ifile}_${datesst}.grib2 ] && [ ! -s ${Dir_ETAwrk_SST}/noaa_sst_${datesst2}.bin ] ; then
   ${Dir_util}/wgrib2 -s -YY  -d 1 -order we:ns ${Dir_SST}/${ifile}_${datesst}.grib2  \
   -ieee ${Dir_ETAwrk_SST}/noaa_sst_${datesst2}.bin
fi
exit
