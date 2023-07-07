#!/bin/ksh -x

Run_Date=$1
Datef=$2
export Dir_scr=DIRROOT
export Dir_home=`dirname ${Dir_scr}`
export Dir_datain=`dirname ${Dir_home}`
export ModelDrive=`basename ${Dir_scr}`
export Dir_wrk=${Dir_datain}/sst/${DirModel}/${Run_Date}
export Dir_ETAwrk_SST=${Dir_datain}/sst/ETAwrk/${DirModel}/${Run_Date}
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
while ((${YY}${MM}<=${Yf}${Mf})); do
  cdo select,name=tos -addc,273.15 ${Dir_wrk}/besm_tos_${YY}_${MM}.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK.nc
  cdo remapbil,besm_sst_grid0.25.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK.nc  ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK_remap.nc
  cdo -sellonlatbox,0,360,90,-90 ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK_remap.nc ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK_0-360_90_-90.nc
  cdo -f grb -copy ${Dir_wrkt}/besm_tos_${YY}_${MM}_TempK_0-360_90_-90.nc ${Dir_wrkt}/besm_sst_${YY}_${MM}_TempK_0-360_90_-90.grib
  Ntimes[$YYMM]=`cdo -ntime ${Dir_wrk}/besm_tos_${YY}_${MM}.nc`
  echo ${Ntimes[$YYMM]}
  let TNTimes=$TNTimes+${Ntimes[$YYMM]}
  let cont=${cont}+32 # soma 32 para definir o proximo mes
  YY=`${Dir_util}/caldate.3.0 ${Run_Date} + ${cont}d "yyyy"`
  MM=`${Dir_util}/caldate.3.0 ${Run_Date} + ${cont}d "mm"`
  YYMM=${YY}${MM}
done

echo ${TNTimes}
echo ${TNTimes} >${Dir_ETAwrk_SST}/NTimes.txt

cont=0
YY=${Yi}
MM=${Mi}
YYMM=${YY}${MM}
#JLG O arquivo ja contem o dado da analise
#${Dir_util}/wgrib -d 1 ${Dir_wrkt}/besm_sst_${YY}_${MM}_TempK_0-360_90_-90.grib -nh -bin -o ${Dir_wrkt}/wfile_0-360_-90_90.bin
i=1
while (($i<=$Ntimes[$YYMM]));do
  ${Dir_util}/wgrib -d $i ${Dir_wrkt}/besm_sst_${YY}_${MM}_TempK_0-360_90_-90.grib -nh -bin -append -o ${Dir_wrkt}/wfile_0-360_-90_90.bin
  let i=$i+1
done

 let cont=${cont}+32
  YY=`${Dir_util}/caldate.3.0 ${Run_Date} + ${cont}d "yyyy"`
  MM=`${Dir_util}/caldate.3.0 ${Run_Date} + ${cont}d "mm"`
  YYMM=${YY}${MM}

while ((${YY}${MM}<=${Yf}${Mf})); do
 ${Dir_util}/wgrib -d 1 ${Dir_wrkt}/besm_sst_${YY}_${MM}_TempK_0-360_90_-90.grib -nh -bin -append -o ${Dir_wrkt}/wfile_0-360_-90_90.bin
 i=2
 while (($i<=${Ntimes[$YYMM]}));do
   ${Dir_util}/wgrib -d $i ${Dir_wrkt}/besm_sst_${YY}_${MM}_TempK_0-360_90_-90.grib -nh -bin -append -o ${Dir_wrkt}/wfile_0-360_-90_90.bin  
   let i=$i+1
 done

 let cont=${cont}+32
  YYS=${YY}
  MMS=${MM}
  YY=`${Dir_util}/caldate.3.0 ${Run_Date} + ${cont}d "yyyy"`
  MM=`${Dir_util}/caldate.3.0 ${Run_Date} + ${cont}d "mm"`
  YYMM=${YY}${MM}
done
  let i=$i-1
  ${Dir_util}/wgrib -d $i ${Dir_wrkt}/besm_sst_${YYS}_${MMS}_TempK_0-360_90_-90.grib -nh -bin -append -o ${Dir_wrkt}/wfile_0-360_-90_90.bin 

  mv ${Dir_wrkt}/wfile_0-360_-90_90.bin ${Dir_ETAwrk_SST}/${DirModel}_sst_${Run_Date}.bin 
  rm -Rf ${Dir_wrkt}
  exit
