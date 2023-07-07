#!/bin/ksh -x
#
#
# Verifica se foi passado algum parametro
hostname
export hh=${1}
if (($#<3)) ; then
export Run_Date=`date "+"%Y%m%d"`${hh}
else
export Run_Date=${3}${hh}
fi

# VARIAVEIS
#export Run_Date=2017103000
export CFS_Member=${2}
export Dir_home=/pesq/grpeta/dsk003/scripts/CFS
export Dir_wrk=/dados/grpeta/dsk003/dados/CFS
export Dir_SST=/dados/grpeta/dsk003/dados/SST
export Dir_scr=${Dir_home}/scripts
export Dir_CFS_out=${Dir_wrk}/data/CFS/output/${Run_Date}_M${CFS_Member}
export Dir_util=${Dir_home}/util
export data_grib=${Dir_CFS}/data/CFS
export Fct=3348
export FInitBC=CFS
export InitBC=6


cd ${Dir_scr}

typeset -Z4 ifct ftime

ifct=000
ihora=`echo ${Run_Date} |cut -c 9-10`

export data_grib=${Dir_wrk}/data/${FInitBC}
export data_sst=${Dir_wrk}/data/grib
export sst_get_date=`${Dir_util}/caldate.3.0 ${Run_Date} - 1d 'yyyymmdd'`
export hh=`echo ${Run_Date} | cut -c9-10`
export date_cfs=`echo ${Run_Date} | cut -c1-8`
export date_sst=`echo ${Run_Date} | cut -c1-8`


echo "CFS DATA  " ${Run_Date}

#ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.2008081300/

datadir=/data/nccf/com/cfs/prod/cfs/cfs.${date_cfs}/${hh}/6hrly_grib_${CFS_Member}
datadirsst=/data/nccf/com/gfs/prod
server=ftpprd.ncep.noaa.gov

fsst=rtgssthr_grb_0.083
ftime=${ifct}

if [ ! -d ${data_grib}/${Run_Date}_M${CFS_Member} ] ; then
  mkdir -p  ${data_grib}/${Run_Date}_M${CFS_Member}
fi


if [[ ! -s ${Dir_SST}/sst_Eta_${date_sst}.bin ]] ; then
   cd ${Dir_SST}
   wget -c http://$server$datadirsst/sst.${sst_get_date}/${fsst}.grib2
   mv -v ${Dir_SST}/${fsst}.grib2 ${Dir_SST}/${fsst}_${date_sst}.grib2
   ${Dir_home}/util/wgrib2 -s -YY  -d 1 -order we:ns ${Dir_SST}/${fsst}_${date_sst}.grib2 -ieee ${Dir_SST}/sst_Eta_${date_sst}.bin
   echo "passou wgrib2"
   chmod 755 ${Dir_SST}/sst_Eta_${date_sst}.bin
fi

cd ${data_grib}/${Run_Date}_M${CFS_Member}
    
    while ((ftime<=${Fct})) ; do
    datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
#for ftime in 00 06 12 18 24 30 36 42 48;do
#  ifile=gfs.t${ihora}z.pgrb2.0p25.f${ftime}.${Run_Date}.grib2
  ifile=pgbf${datefct}.${CFS_Member}.${Run_Date}.grb2
  ifilesolo=flxf${datefct}.${CFS_Member}.${Run_Date}.grb2
  
  if [[ ! -s ${data_grib}/${Run_Date}/${ifile} ]] && [[ ! -s ${data_sst}/${fsst}_${Run_Date}.grib2 ]] ; then
    echo $server$datadir/${ifile} $server$datadir/${ifilesolo}
    echo "http://$server$datadir/${ifile}" >> ${data_grib}/${Run_Date}_M${CFS_Member}/wget_${CFS_Member}.list
    echo "http://$server$datadir/${ifilesolo}" >> ${data_grib}/${Run_Date}_M${CFS_Member}/wget_${CFS_Member}.list
   fi
        let ftime=$ftime+$InitBC

    done
 
   cat ${data_grib}/${Run_Date}_M${CFS_Member}/wget_${CFS_Member}.list | xargs -n 1 -P 16 wget -q -c
   
   ${Dir_scr}/CFS_deco.ksh ${hh} ${CFS_Member} ${date_cfs}

   cd ${data_grib}/${Run_Date}_M${CFS_Member}
   tar -zcvf ${Run_Date}_M${CFS_Member}.tgz pgbf*${CFS_Member}.${Run_Date}.grb2 flxf*${CFS_Member}.${Run_Date}.grb2
 
#  scp -r ${Dir_CFS_out} jorge.gomes@tupa:/scratchout/grupos/eta/home/jorge.gomes/CFS/.
#  ssh  jorge.gomes@tupa 'chmod -R 770 /scratchout/grupos/eta/home/jorge.gomes/CFS/'${Run_Date}'_M'${CFS_Member}
#  scp  ${Dir_SST}/sst_Eta_${date_sst}.bin jorge.gomes@tupa:/scratchout/grupos/eta/home/jorge.gomes/SST/.
#  ssh  jorge.gomes@tupa 'chmod -R 770 /scratchout/grupos/eta/home/jorge.gomes/SST/sst_Eta_${date_sst}.bin'
   
   
   
   
