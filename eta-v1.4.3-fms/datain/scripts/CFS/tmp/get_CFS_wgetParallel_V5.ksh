#!/bin/ksh -x
#
#
# Verifica se foi passado algum parametro
hostname
export hh=${1}
if (($#<4)) ; then
export Run_Date=`date "+"%Y%m%d"`${hh}
else
export Run_Date=${5}${hh}
fi

# VARIAVEIS
export FInitBC=CFS
export CFS_Member=${2}
export Dir_home=/dados/grpeta/dsk003/scripts/CFS
export Dir_wrk=/dados/grpeta/dsk003/dados/CFS
export Dir_SST=/dados/grpeta/dsk003/dados/SST
export Dir_scr=${Dir_home}/scripts
export Dir_CFS_out=${Dir_wrk}/output/${Run_Date}_M${CFS_Member}
export Dir_util=${Dir_home}/util
export Fct=${4}
export InitBC=6
export IntSST=24
export ifct=${3}

cd ${Dir_scr}

typeset -Z4 ifct ftime Fct

#ifct=000
ihora=`echo ${Run_Date} |cut -c 9-10`

export data_grib=${Dir_wrk}

export yyyy=`echo ${Run_Date} | cut -c1-4`
export mm=`echo ${Run_Date} | cut -c5-6`
export dd=`echo ${Run_Date} | cut -c7-8`
export hh=`echo ${Run_Date} | cut -c9-10`
export date_cfs=`echo ${Run_Date} | cut -c1-8`
export date_sst=`echo ${Run_Date} | cut -c1-8`


echo "CFS DATA  " ${Run_Date}


#https://www.ncei.noaa.gov/data/climate-forecast-system/access/operational-9-month-forecast/6-hourly-by-pressure/2021/202109/20210909/2021090900/
dir_pgbf=/6-hourly-by-pressure/${yyyy}/${yyyy}${mm}/${yyyy}${mm}${dd}/${yyyy}${mm}${dd}${hh}
dir_flxf=/6-hourly-flux/${yyyy}/${yyyy}${mm}/${yyyy}${mm}${dd}/${yyyy}${mm}${dd}${hh}
dir_ocnf=/6-hourly-ocean/${yyyy}/${yyyy}${mm}/${yyyy}${mm}${dd}/${yyyy}${mm}${dd}${hh}


datadir=/data/nccf/com/cfs/prod/cfs/cfs.${date_cfs}/${hh}/6hrly_grib_${CFS_Member}
datadirsst=/data/nccf/com/gfs/prod
server=www.ncei.noaa.gov/data/climate-forecast-system/access/operational-9-month-forecast

ftime=${ifct}

if [ ! -d ${data_grib}/${Run_Date}_M${CFS_Member} ] ; then
  mkdir -p  ${data_grib}/${Run_Date}_M${CFS_Member}
else
  rm -f ${data_grib}/${Run_Date}_M${CFS_Member}/wget_pgbf_${CFS_Member}.list 
  rm -f ${data_grib}/${Run_Date}_M${CFS_Member}/pgbf_*
  rm -f ${data_grib}/${Run_Date}_M${CFS_Member}/wget_flxf_${CFS_Member}.list
  rm -f ${data_grib}/${Run_Date}_M${CFS_Member}/flxf_*
  rm -f ${data_grib}/${Run_Date}_M${CFS_Member}/wget_ocnf_${CFS_Member}.list 
  rm -f ${data_grib}/${Run_Date}_M${CFS_Member}/ocnf_*
fi


cd ${data_grib}/${Run_Date}_M${CFS_Member}
    
    while ((ftime<=${Fct})) ; do
    datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
#for ftime in 00 06 12 18 24 30 36 42 48;do
#  ifile=gfs.t${ihora}z.pgrb2.0p25.f${ftime}.${Run_Date}.grib2
  ifile=pgbf${datefct}.${CFS_Member}.${Run_Date}.grb2
  ifilesolo=flxf${datefct}.${CFS_Member}.${Run_Date}.grb2
  
  if [[ ! -s ${data_grib}/${Run_Date}/${ifile} ]] ; then
    echo $server$dir_pgbf/${ifile} $server$dir_flxf/${ifilesolo}
    echo "https://$server$dir_pgbf/${ifile}" >> ${data_grib}/${Run_Date}_M${CFS_Member}/wget_pgbf_${CFS_Member}.list
    echo "https://$server$dir_flxf/${ifilesolo}" >> ${data_grib}/${Run_Date}_M${CFS_Member}/wget_flxf_${CFS_Member}.list
   fi
        let ftime=$ftime+$InitBC

    done
    ftime=24   
    while ((ftime<=${Fct})) ; do
      datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
      ifilesst=ocnf${datefct}.${CFS_Member}.${Run_Date}.grb2
      if [[ ! -s ${data_grib}/${Run_Date}/${ifilesst} ]] ; then
        echo "https://$server$dir_ocnf/${ifilesst}" >> ${data_grib}/${Run_Date}_M${CFS_Member}/wget_ocnf_${CFS_Member}.list
      fi
        let ftime=$ftime+$IntSST
    done
 
cd ${data_grib}/${Run_Date}_M${CFS_Member}
split  -l20 ${data_grib}/${Run_Date}_M${CFS_Member}/wget_pgbf_${CFS_Member}.list pgbf_
for fl in $(ls  pgbf_*)
do
   print "Current file : $fl"
   cat ${data_grib}/${Run_Date}_M${CFS_Member}/${fl}  | xargs -n 1 -P 10 wget -q -c -t 3
done

split  -l20 ${data_grib}/${Run_Date}_M${CFS_Member}/wget_flxf_${CFS_Member}.list flxf_
for fl in $(ls  flxf_*)
do
   print "Current file : $fl"
   cat ${data_grib}/${Run_Date}_M${CFS_Member}/${fl}  | xargs -n 1 -P 10 wget -q -c -t 3 
done

split  -l20 ${data_grib}/${Run_Date}_M${CFS_Member}/wget_ocnf_${CFS_Member}.list ocnf_
for fl in $(ls  ocnf_*)
do
   print "Current file : $fl"
   cat ${data_grib}/${Run_Date}_M${CFS_Member}/${fl}  | xargs -n 1 -P 10 wget -q -c -t 3
done
 
${Dir_scr}/CFS_deco_V5.ksh ${hh} ${CFS_Member} ${ifct} ${Fct} ${date_cfs}
   
   
