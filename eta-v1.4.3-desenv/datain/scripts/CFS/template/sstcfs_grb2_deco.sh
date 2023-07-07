#!/bin/ksh -x
#
# Verifica se foi passado algum parametro
if (($#<5)) ; then
export Run_Date=`date "+%Y%m%d"`$1
else
export Run_Date=${5}${1}
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
export Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}.${CFS_Member}/${Run_Date}
export Dir_SST_ETAwrk=${Dir_datain}/sst/ETAwrk/${ModelDrive}.${CFS_Member}/${Run_Date}
export Dir_util=${Dir_datain}/util
export Dir_dprep_exe=${Dir_datain}/dprep/exe
export InitBC=6
export IntSST=24
export FInitBC=${ModelDrive}.${CFS_Member}

mkdir -p ${Dir_SST_ETAwrk}
me=$(whoami)
echo "I am $me."
ofilesst=${FInitBC}_sst_${Run_Date}.bin
echo 'remove old files!!!'
rm -f  ${Dir_SST_ETAwrk}/${ofilesst}
touch  ${Dir_SST_ETAwrk}/${ofilesst}
echo 'OPEN ' ${Dir_SST_ETAwrk}/${ofilesst}'  FILE'

datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${Fcti}hr 'yyyymmddhh'`
ifilesst=ocnf${datefct}.${CFS_Member}.${Run_Date}.grb2

line=$(${Dir_util}/wgrib2 -ctl_inv ${Dir_wrk}'/'${ifilesst} | grep "POT 160,5 0,0,2")

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -s -YY  -d $nrec -order we:ns ${Dir_wrk}'/'${ifilesst} -ncpu 4 -append -ieee ${Dir_SST_ETAwrk}/${ofilesst} >> ${Dir_SST_ETAwrk}/log
typeset -Z4 NTTimes
NTTimes=1
ftime=${Fcti}
while [ ${ftime} -le ${Fct} ];do
      datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${ftime}hr 'yyyymmddhh'`
      ifilesst=ocnf${datefct}.${CFS_Member}.${Run_Date}.grb2

echo '==================================================='
echo 'file in process '${Dir_wrk}'/'${ifilesst}
echo '==================================================='


line=$(${Dir_util}/wgrib2 -ctl_inv ${Dir_wrk}'/'${ifilesst} | grep "POT 160,5 0,0,2")

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -s -YY  -d $nrec -order we:ns ${Dir_wrk}'/'${ifilesst} -ncpu 4 -append -ieee ${Dir_SST_ETAwrk}/${ofilesst} >> ${Dir_SST_ETAwrk}/log


echo 'DONE WITH ' ${Dir_SST_ETAwrk}'/'${ifilesst}'  FILE'

ftime=`expr $ftime + ${IntSST} `
done # while [ $dhr -le $dend ]
echo "${NTTimes}" > ${Dir_SST_ETAwrk}/NSSTFLDS.txt

exit
