#!/bin/ksh -x
#
#
# Verifica se foi passado algum parametro
hostname
export hh=${1}
if (($#<2)) ; then
export Run_Date=`date "+"%Y%m%d"`${hh}
else
export Run_Date=${2}${hh}
fi

# VARIAVEIS
#export Run_Date=2017103000
export Dir_home=/pesq/grpeta/dsk003/scripts/CFS
export Dir_wrk=/dados/grpeta/dsk003/dados/CFS
export Dir_SST=/dados/grpeta/dsk003/dados/SST
export Dir_scr=${Dir_home}/scripts
export Dir_util=${Dir_home}/util


cd ${Dir_scr}


ifct=000
ihora=`echo ${Run_Date} |cut -c 9-10`

export sst_get_date=`${Dir_util}/caldate.3.0 ${Run_Date} - 1d 'yyyymmdd'`
export hh=`echo ${Run_Date} | cut -c9-10`
export date_sst=`echo ${Run_Date} | cut -c1-8`


echo "DATA  " ${Run_Date}

#ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.2008081300/

datadirsst=/data/nccf/com/gfs/prod
server=ftpprd.ncep.noaa.gov

fsst=rtgssthr_grb_0.083


if [[ ! -s ${Dir_SST}/sst_Eta_${date_sst}.bin ]] ; then
   cd ${Dir_SST}
   wget -c http://$server$datadirsst/sst.${sst_get_date}/${fsst}.grib2
   mv -v ${Dir_SST}/${fsst}.grib2 ${Dir_SST}/${fsst}_${date_sst}.grib2
   ${Dir_home}/util/wgrib2 -s -YY  -d 1 -order we:ns ${Dir_SST}/${fsst}_${date_sst}.grib2 -ieee ${Dir_SST}/sst_Eta_${date_sst}.bin
   echo "passou wgrib2"
   chmod 755 ${Dir_SST}/sst_Eta_${date_sst}.bin
else
   echo "SST ${Run_Date} available"   
fi   
   
   
   
