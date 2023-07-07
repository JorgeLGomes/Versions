#!/bin/ksh -x

Run_Date=$1
Datef=$2
export Dir_scr=DIRROOT
export Dir_home=`dirname ${Dir_scr}`
export Dir_datain=`dirname ${Dir_home}`
export ModelDrive=`basename ${Dir_scr}`
export Dir_wrk=${Dir_datain}/sst/${ModelDrive}/${Run_Date}
export Dir_ETAwrk_SST=${Dir_datain}/sst/ETAwrk/${ModelDrive}.mean/${Run_Date}
export Dir_wrkt=${Dir_ETAwrk_SST}/WRK
export Dir_util=${Dir_datain}/util

if [ ! -d ${Dir_ETAwrk_SST} ] ; then
  mkdir -p  ${Dir_ETAwrk_SST}
fi
if [ ! -d ${Dir_wrkt} ] ; then
  mkdir -p  ${Dir_wrkt}
fi
Yi=`echo ${Run_Date}| cut -c1-4`
Mi=`echo ${Run_Date}| cut -c5-6`

Yf=`echo ${Datef}| cut -c1-4`
Mf=`echo ${Datef}| cut -c5-6`

typeset -Z2 Mi Mf
typeset -Z4 TNTimes

TNTimes=1
cont=0
YY=${Yi}
MM=${Mi}
YYMM=${YY}${MM}
cdo select,name=tos -addc,273.15 ${Dir_wrk}/besm_tos_${YY}_${MM}.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK.nc
cdo chname,tos,sst ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_CVname.nc
cdo remapbil,besm_sst_grid0.25.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_CVname.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_remap.nc
cdo -sellonlatbox,0,360,90,-90 ${Dir_wrkt}/besm_tos_${YY}_${MM}_remap.nc ${Dir_wrkt}/besm_tos_${YY}_${MM}_0-360_90_-90.nc
cdo monmean ${Dir_wrkt}/besm_tos_${YY}_${MM}_0-360_90_-90.nc ${Dir_wrkt}/besm_tos_${YY}_${MM}_MonMean.nc
cdo -f grb -copy ${Dir_wrkt}/besm_tos_${YY}_${MM}_MonMean.nc ${Dir_wrkt}/besm_sst_${YY}_${MM}.grib
${Dir_util}/wgrib -d 1 ${Dir_wrkt}/besm_sst_${YY}_${MM}.grib -nh -bin -append -o ${Dir_wrkt}/wfile.bin
while ((${YY}${MM}<=${Yf}${Mf})); do
  cdo select,name=tos -addc,273.15 ${Dir_wrk}/besm_tos_${YY}_${MM}.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK.nc
  cdo chname,tos,sst ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_CVname.nc 
  cdo remapbil,besm_sst_grid0.25.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_CVname.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_remap.nc
  cdo -sellonlatbox,0,360,90,-90 ${Dir_wrkt}/besm_tos_${YY}_${MM}_remap.nc ${Dir_wrkt}/besm_tos_${YY}_${MM}_0-360_90_-90.nc
  cdo monmean ${Dir_wrkt}/besm_tos_${YY}_${MM}_0-360_90_-90.nc ${Dir_wrkt}/besm_tos_${YY}_${MM}_MonMean.nc
  cdo -f grb -copy ${Dir_wrkt}/besm_tos_${YY}_${MM}_MonMean.nc ${Dir_wrkt}/besm_sst_${YY}_${MM}.grib
  ${Dir_util}/wgrib -d 1 ${Dir_wrkt}/besm_sst_${YY}_${MM}.grib -nh -bin -append -o ${Dir_wrkt}/wfile.bin
  let cont=${cont}+32 # soma 32 para definir o proximo mes
  YYS=${YY}
  MMS=${MM}
  YY=`${Dir_util}/caldate.3.0 ${Run_Date} + ${cont}d "yyyy"`
  MM=`${Dir_util}/caldate.3.0 ${Run_Date} + ${cont}d "mm"`
  YYMM=${YY}${MM}
  let TNTimes=${TNTimes}+1
done

  let i=$i-1
  ${Dir_util}/wgrib -d 1 ${Dir_wrkt}/besm_sst_${YYS}_${MMS}.grib -nh -bin -append -o ${Dir_wrkt}/wfile.bin 
  let TNTimes=${TNTimes}+1
  echo ${TNTimes} >${Dir_ETAwrk_SST}/NTimes.txt
  mv ${Dir_wrkt}/wfile.bin ${Dir_ETAwrk_SST}/${ModelDrive}.mean_sst_${Run_Date}.bin 
#  rm -Rf ${Dir_wrkt}
exit 0
