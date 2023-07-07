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
Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
Dir_SST_ETAwrk=${Dir_datain}/sst/ETAwrk/noaa/${Run_Date}
Dir_util=${Dir_datain}/util
InitBC=6
FInitBC=${ModelDrive}

if [ ! -d ${Dir_SST_ETAwrk} ] ; then
  mkdir -p ${Dir_SST_ETAwrk}
fi

if [ ! -d ${Dir_ETAwrk} ] ; then
  mkdir -p ${Dir_ETAwrk}
  mkdir -p ${Dir_SST_ETAwrk}
else
  rm -f ${Dir_ETAwrk}/wget.list 
fi

cd ${Dir_scr}

FctiF=`printf "%03d" "${Fcti}"`
FctF=`printf "%03d" "${Fct}"`

export yyyy=`echo ${Run_Date} | cut -c1-4`
export mm=`echo ${Run_Date} | cut -c5-6`
export dd=`echo ${Run_Date} | cut -c7-8`
export hh=`echo ${Run_Date} | cut -c9-10`

echo "GFS DATA  " ${Run_Date}

server=http://ftp1.cptec.inpe.br
datadir=pesquisa/grpeta/VII-WorkEta/InputData/Weather
datadir2=sst/ETAwrk/noaa/${Run_Date}
ifile1=noaa_sst_${Run_Date}.bin.xz
cd ${Dir_SST_ETAwrk}
wget -c $server/${datadir}/${datadir2}/${ifile1}
xz -f -d ${ifile1}

datadir2=atmos/ETAwrk/${ModelDrive}/${Run_Date}
cd ${Dir_ETAwrk}
    
ftime=${Fcti}
ftimeF=`printf "%06d" "${iftime}"`
while ((ftime<=${Fct})) ; do
  datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
  ifile1=${ModelDrive}_${Run_Date}.${ftimeF}.ETAwrk.xz
  ifile2=${ModelDrive}_${Run_Date}.${ftimeF}.ETAwrk.gdsinfo.xz
  echo $server/${datadir}/${datadir2}/${ifile1}
  echo "$server/${datadir}/${datadir2}/${ifile1}" >> ${Dir_ETAwrk}/wget.list
  echo "$server/${datadir}/${datadir2}/${ifile2}" >> ${Dir_ETAwrk}/wget.list
  echo "${ifile1}" >> ${Dir_ETAwrk}/xz.list
  echo "${ifile2}" >> ${Dir_ETAwrk}/xz.list
  let ftime=$ftime+$InitBC
  ftimeF=`printf "%06d" "${ftime}"`
done

cd ${Dir_ETAwrk}
print "Current file : wget.list"
cat ${Dir_ETAwrk}/wget.list  | xargs -n 1 -P 10 wget -q -c -t 3
cat ${Dir_ETAwrk}/xz.list    | xargs -n 1 -P 10 xz -f -T 0 -d 
exit
