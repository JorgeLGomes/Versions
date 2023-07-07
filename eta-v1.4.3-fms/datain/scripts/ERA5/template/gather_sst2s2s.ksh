#!/bin/ksh -x

Datei=$1
Datef=$2
InitBC=$3
export Dir_scr=DIRROOT
export Dir_home=`dirname ${Dir_scr}`
export Dir_datain=`dirname ${Dir_home}`
export ModelDrive=`basename ${Dir_scr}`
export Dir_ETAwrk_SST=${Dir_datain}/sst/ETAwrk/${ModelDrive}
export Dir_util=${Dir_datain}/util
export Dir_ETAwrk_SST_S2S=${Dir_datain}/sst/ETAwrk/${ModelDrive}_S2S/${Datei}

mkdir ${Dir_ETAwrk_SST_S2S}
typeset -Z2 Mi Mf
typeset -Z4 TNTimes
Run_Date=${Datei}
cat ${Dir_ETAwrk_SST}/${Run_Date}/ERA5_sst_${Run_Date}.bin > ${Dir_ETAwrk_SST_S2S}/ERA5_S2S_sst_${Datei}.bin
Run_Date=`${Dir_util}/caldate.3.0 ${Run_Date} + ${InitBC}hr 'yyyymmddhh'`
TNTimes=1
cont=0
while [ "${Run_Date}" != "${Datef}" ] ; do
  cat ${Dir_ETAwrk_SST}/${Run_Date}/ERA5_sst_${Run_Date}.bin >> ${Dir_ETAwrk_SST_S2S}/ERA5_S2S_sst_${Datei}.bin
  Run_Date=`${Dir_util}/caldate.3.0 ${Run_Date} + ${InitBC}hr 'yyyymmddhh'`
  let TNTimes=${TNTimes}+1
done
echo ${TNTimes}
echo ${TNTimes} >${Dir_ETAwrk_SST_S2S}/NTimes.txt
exit
