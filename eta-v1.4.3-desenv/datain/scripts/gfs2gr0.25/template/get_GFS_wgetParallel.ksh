#!/bin/ksh -x
#
#
# Verifica se foi passado algum parametro
hostname
export hh=${1}
if (($#<4)) ; then
export Run_Date=`date "+%Y%m%d"`${hh}
else
export Run_Date=${4}${hh}
fi

export Fcti=${2}
export Fct=${3}

# VARIAVEIS
Dir_scr=DIRROOT
Dir_home=`dirname ${Dir_scr}`
Dir_datain=`dirname ${Dir_home}`
ModelDrive=`basename ${Dir_scr}`
Dir_wrk=${Dir_datain}/atmos/${ModelDrive}/${Run_Date}
Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
Dir_util=${Dir_datain}/util
InitBC=6
FInitBC=${ModelDrive}

if [ ! -d ${Dir_wrk} ] ; then
  mkdir -p  ${Dir_wrk}
  mkdir -p ${Dir_ETAwrk}
else
  rm -f ${Dir_wrk}/wget_pgrb2.list 
fi

cd ${Dir_scr}

FctiF=`printf "%03d" "${Fcti}"`
FctF=`printf "%03d" "${Fct}"`

export yyyy=`echo ${Run_Date} | cut -c1-4`
export mm=`echo ${Run_Date} | cut -c5-6`
export dd=`echo ${Run_Date} | cut -c7-8`
export hh=`echo ${Run_Date} | cut -c9-10`

echo "GFS DATA  " ${Run_Date}

#https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20220112/00/atmos/
datadir=pub/data/nccf/com/gfs/prod
server=https://nomads.ncep.noaa.gov
dir_pgrb2=gfs.${yyyy}${mm}${dd}/${hh}/atmos
fpref=gfs.t${hh}z.pgrb2.0p25.f     

cd ${Dir_wrk}
    
ftime=${Fcti}
ftimeF=`printf "%03d" "${iftime}"`
while ((ftime<=${Fct})) ; do
  datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
  ifile=gfs.t${hh}z.pgrb2.0p25.f${ftimeF}
  if [[ ! -s ${Dir_wrk}/${ifile} ]] ; then
    echo $server/${datadir}/${dir_pgrb2}/${ifile}
    echo "$server/${datadir}/${dir_pgrb2}/${ifile}" >> ${Dir_wrk}/wget_pgrb2.list
  fi
  let ftime=$ftime+$InitBC
  ftimeF=`printf "%03d" "${ftime}"`
done

cd ${Dir_wrk}
print "Current file : wget_pgrb2.list"
cat ${Dir_wrk}/wget_pgrb2.list  | xargs -n 1 -P 10 wget --no-check-certificate -q -c -t 3
exit
