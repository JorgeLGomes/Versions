#!/bin/ksh -x
#
#
# Verifica se foi passado algum parametro
hostname
export hh=${1}
if (($#<5)) ; then
export Run_Date=`date "+%Y%m%d"`${hh}
else
export Run_Date=${5}${hh}
fi

export CFS_Member=${2}
export Fcti=${3} 
export Fct=${4}

# VARIAVEIS
export Dir_scr=DIRROOT
export Dir_home=`dirname ${Dir_scr}`
export Dir_datain=`dirname ${Dir_home}`
export ModelDrive=`basename  ${Dir_scr}`
export Dir_wrk=${Dir_datain}/atmos/${ModelDrive}.${CFS_Member}/${Run_Date}
export Dir_util=${Dir_datain}/util
export InitBC=6
export IntSST=24

cd ${Dir_scr}

typeset -Z4 Fcti ftime Fct

ihora=`echo ${Run_Date} |cut -c 9-10`

export yyyy=`echo ${Run_Date} | cut -c1-4`
export mm=`echo ${Run_Date} | cut -c5-6`
export dd=`echo ${Run_Date} | cut -c7-8`
export hh=`echo ${Run_Date} | cut -c9-10`
export date_cfs=`echo ${Run_Date} | cut -c1-8`
export date_sst=`echo ${Run_Date} | cut -c1-8`


echo "${ModelDrive} DATA  " ${Run_Date}


#https://www.ncei.noaa.gov/data/climate-forecast-system/access/operational-9-month-forecast/6-hourly-by-pressure/2021/202109/20210909/2021090900/
dir_pgbf=/6-hourly-by-pressure/${yyyy}/${yyyy}${mm}/${yyyy}${mm}${dd}/${yyyy}${mm}${dd}${hh}
dir_flxf=/6-hourly-flux/${yyyy}/${yyyy}${mm}/${yyyy}${mm}${dd}/${yyyy}${mm}${dd}${hh}
dir_ocnf=/6-hourly-ocean/${yyyy}/${yyyy}${mm}/${yyyy}${mm}${dd}/${yyyy}${mm}${dd}${hh}


datadir=/data/nccf/com/cfs/prod/cfs/cfs.${date_cfs}/${hh}/6hrly_grib_${CFS_Member}
datadirsst=/data/nccf/com/gfs/prod
server=www.ncei.noaa.gov/data/climate-forecast-system/access/operational-9-month-forecast

ftime=${Fcti}

if [ -d ${Dir_wrk} ] ; then
  rm -f ${Dir_wrk}/wget_pgbf_${CFS_Member}.list 
  rm -f ${Dir_wrk}/pgbf_*
  rm -f ${Dir_wrk}/wget_flxf_${CFS_Member}.list
  rm -f ${Dir_wrk}/flxf_*
  rm -f ${Dir_wrk}/wget_ocnf_${CFS_Member}.list 
  rm -f ${Dir_wrk}/ocnf_*
fi


cd ${Dir_wrk}
    
while ((ftime<=${Fct})) ; do
  datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
  ifile=pgbf${datefct}.${CFS_Member}.${Run_Date}.grb2
  ifilesolo=flxf${datefct}.${CFS_Member}.${Run_Date}.grb2
  
  if [[ ! -s ${Dir_wrk}/${Run_Date}/${ifile} ]] ; then
    echo $server$dir_pgbf/${ifile} $server$dir_flxf/${ifilesolo}
    echo "https://$server$dir_pgbf/${ifile}" >> ${Dir_wrk}/wget_pgbf_${CFS_Member}.list
    echo "https://$server$dir_flxf/${ifilesolo}" >> ${Dir_wrk}/wget_flxf_${CFS_Member}.list
   fi
        let ftime=$ftime+$InitBC

    done
    ftime=24   
    while ((ftime<=${Fct})) ; do
      datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
      ifilesst=ocnf${datefct}.${CFS_Member}.${Run_Date}.grb2
      if [[ ! -s ${Dir_wrk}/${ifilesst} ]] ; then
        echo "https://$server$dir_ocnf/${ifilesst}" >> ${Dir_wrk}/wget_ocnf_${CFS_Member}.list
      fi
        let ftime=$ftime+$IntSST
    done
 
cd ${Dir_wrk}
split  -l20 ${Dir_wrk}/wget_pgbf_${CFS_Member}.list pgbf_
for fl in $(ls  pgbf_*)
do
   print "Current file : $fl"
   cat ${Dir_wrk}/${fl}  | xargs -n 1 -P 10 wget --no-check-certificate -q -c -t 3
done

split  -l20 ${Dir_wrk}/wget_flxf_${CFS_Member}.list flxf_
for fl in $(ls  flxf_*)
do
   print "Current file : $fl"
   cat ${Dir_wrk}/${fl}  | xargs -n 1 -P 10 wget --no-check-certificate -q -c -t 3 
done

split  -l20 ${Dir_wrk}/wget_ocnf_${CFS_Member}.list ocnf_
for fl in $(ls  ocnf_*)
do
   print "Current file : $fl"
   cat ${Dir_wrk}/${fl}  | xargs -n 1 -P 10 wget --no-check-certificate -q -c -t 3
done
 
#${Dir_scr}/CFS_deco_V5.ksh ${hh} ${CFS_Member} ${ifct} ${Fct} ${date_cfs}
   
   
