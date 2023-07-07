#!/bin/ksh -x
#
# Verifica se foi passado algum parametro
if (($#<3)) ; then
export Run_Date=`date "+"%Y%m%d"`$1
else
export Run_Date=${3}${1}
fi

# VARIAVEIS
#export Run_Date=2017103000
export CFS_Member=${2}
export Dir_home=/dados/grpeta/dsk003/scripts/CFS
export Dir_wrk=/dados/grpeta/dsk003/dados/CFS
export Dir_SST=/dados/grpeta/dsk003/dados/SST
export Dir_scr=${Dir_home}/scripts
export Dir_wrk_out=${Dir_wrk}/data/CFS/output/${Run_Date}_M${CFS_Member}
#export Dir_util=${Dir_home}/util
export Dir_util=/scratchin/grupos/grpeta/outros/unified/Eta_support_data/util
export data_grib=${Dir_wrk}/data/CFS
export Fct=1440
export FInitBC=CFS
export Int=24

export GRBDIR=${data_grib}/${Run_Date}_M${CFS_Member}
export WOUTDIR=${Dir_wrk_out}
mkdir -p ${Dir_wrk_out}
me=$(whoami)
echo "I am $me."

rm -f ${WOUTDIR}/log

typeset -Z4 dhr
typeset -Z4 dhri Fct

ifct=024
dhri=${ifct}
dhr=${dhri}

datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${dhri}hr 'yyyymmddhh'`

rm -f  ${WOUTDIR}/SST_${adate}
touch  ${WOUTDIR}/SST_${adate}
echo 'OPEN '${WOUTDIR}'/SST_'${adate}'  FILE'

line=$(${Dir_util}/wgrib2 -v $GRBDIR'/ocnf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "POT" | grep ":5 m")
adate=`echo ${line}|cut -d"=" -f2|cut -c1-10`
echo $adate > ${WOUTDIR}/adate.txt
echo $adate $Run_Date


if [ $adate -eq $Run_Date ] ; then
  echo "The cfs date is $adate"
else
  echo "Data for $Run_Date isn't available"
  exit 99 
fi

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}


${Dir_util}/wgrib2 -s -YY -d $nrec -order we:ns $GRBDIR'/ocnf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' -ieee  ${WOUTDIR}/SST_${adate} >> ${WOUTDIR}/log


ifct=024

dhri=${ifct}


while [ $dhri -le ${Fct} ];do
echo "passei"

dhr=${dhri}
datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${dhri}hr 'yyyymmddhh'`

echo '==================================================='
echo 'file in process '$GRBDIR'/ocnf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2'
echo '==================================================='


line=$(${Dir_util}/wgrib2 -v $GRBDIR'/ocnf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "POT" | grep ":5 m")
tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}


${Dir_util}/wgrib2 -s -YY -d $nrec -order we:ns $GRBDIR'/ocnf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' -append -ieee ${WOUTDIR}/SST_${adate} >> ${WOUTDIR}/log
dhri=`expr $dhri + ${Int} `

done # while [ $dhr -le $dend ]


exit
