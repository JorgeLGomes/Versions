#!/bin/ksh -x
#
#

# VARIAVEIS
Eta_scr=/lustre_xc50/jorge_gomes/sub-sazonal/sst_noaa/scripts
Eta_home1=/lustre_xc50/jorge_gomes/sub-sazonal/worketa
dataout=/lustre_xc50/jorge_gomes/sub-sazonal/sst_noaa/dataout
cd ${Eta_scr}
Run_Date=$1
Ndays=$2
fsst=rtgssthr_grb_0.083.grib2 
cd ${Eta_run}

ContDays=0
rm ${dataout}/sstgrb.${Run_Date}_${Ndays}
touch ${dataout}/sstgrb.${Run_Date}_${Ndays}
while  [ ${ContDays} -le ${Ndays} ] ; do
  sst_rundate=`${Eta_home1}/util/caldate.3.0 ${Run_Date} + ${ContDays}d 'yyyymmddhh'`
  sst_date=`${Eta_home1}/util/caldate.3.0 ${Run_Date} + ${ContDays}d 'yyyymmdd'`
  dir_sst=/lustre_xc50/ioper/data/external/${sst_rundate}/dataout/NCEP
  if  [[ -s ${dir_sst}/${fsst}.${sst_date} ]] ; then
     ${Eta_home1}/util/wgrib2 -s -YY  -d 1 -order we:ns ${dir_sst}/${fsst}.${sst_date} -ieee ${dataout}/sstgrb.${sst_date}
     cat ${dataout}/sstgrb.${Run_Date}_${Ndays} ${dataout}/sstgrb.${sst_date} >> ${dataout}/sstgrb.${Run_Date}_${Ndays}
  else
     echo "${dir_sst}/${fsst}.${sst_date} file not found!"
     exit 0
  fi
  let ContDays=${ContDays}+1
done

exit 1
