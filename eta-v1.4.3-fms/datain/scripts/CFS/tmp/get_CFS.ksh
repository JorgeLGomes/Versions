#!/bin/ksh -x
#
#
# Verifica se foi passado algum parametro
if (($#<2)) ; then
export Run_Date=`date "+"%Y%m%d"`$1
else
export Run_Date=${2}${1}
fi

# VARIAVEIS
#export Run_Date=2017103000
export CFS_Member=01
export Eta_home=/home/sondamod/CFS
export Eta_wrk=/rede/etap01/sonda/CFS
export Eta_scr=${Eta_home}/scripts
export Eta_util=${Eta_home}/util
export data_grib=${Eta_wrk}/data/CFS
export Fct=264
export FInitBC=CFS
export InitBC=6


cd ${Eta_scr}

typeset -Z3 ifct ftime

ifct=000
ihora=`echo ${Run_Date} |cut -c 9-10`

data_grib=${Eta_wrk}/data/${FInitBC}
data_sst=${Eta_wrk}/data/grib
sst_data=`${Eta_util}/caldate.3.0 ${Run_Date} - 1d 'yyyymmdd'`
hh=`echo ${Run_Date} | cut -c9-10`
date_cfs=`echo ${Run_Date} | cut -c1-8`



echo "CFS DATA  " ${Run_Date}

#ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.2008081300/

#datadir=/data/nccf/com/gfs/prod
datadir=/data/nccf/com/cfs/prod/cfs/cfs.${date_cfs}/${hh}/6hrly_grib_${CFS_Member}
datadirsst=/data/nccf/com/gfs/prod
server=ftpprd.ncep.noaa.gov
#fpref=gfs.t${ihora}z.pgrb2.0p25.f
fsst=rtgssthr_grb_0.083
ftime=${ifct}
if [ ! -d ${data_grib}/${Run_Date} ] ; then
  mkdir -p  ${data_grib}/${Run_Date}
fi
cd ${data_grib}
    while ((ftime<=${Fct})) ; do
    datefct=`${Eta_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
#for ftime in 00 06 12 18 24 30 36 42 48;do
#  ifile=gfs.t${ihora}z.pgrb2.0p25.f${ftime}.${Run_Date}.grib2
  ifile=pgbf${datefct}.${CFS_Member}.${Run_Date}.grb2
  ifilesolo=flxf${datefct}.${CFS_Member}.${Run_Date}.grb2
  
  if [[ ! -s ${data_grib}/${Run_Date}/${ifile} ]] && [[ ! -s ${data_sst}/${fsst}_${Run_Date}.grib2 ]] ; then
       if [ "${ftime}" == "000" ] ; then  
         cd ${data_sst}
         wget -c http://$server$datadirsst/sst.${sst_data}/${fsst}.grib2
         data=`echo ${Run_Date} | cut -c1-8 `
         mv -v ${fsst}.grib2 ${data_sst}/${fsst}_${data}.grib2
         ${Eta_home}/util/wgrib2 -s -YY  -d 1 -order we:ns ${data_sst}/${fsst}_${data}.grib2 -ieee ${data_sst}/sst_Eta_${data}.bin
         echo "passou wgrib2"
         chmod 755 ${data_sst}/sst_Eta_${sst_data}.bin
         cd ${Eta_wrk}
        fi
    echo $server$datadir/${ifile} $server$datadir/${ifilesolo}
    wget -c http://$server$datadir/${ifile}
    wget -c http://$server$datadir/${ifilesolo}

    mv -v  ${ifile} ${data_grib}/${Run_Date}/${ifile}
    mv -v  ${ifilesolo} ${data_grib}/${Run_Date}/${ifilesolo}
    chmod 755 ${data_grib}/${Run_Date}/${ifile}
    chmod 755 ${data_grib}/${Run_Date}/${ifilesolo}
   fi

        let ftime=$ftime+$InitBC

    done

