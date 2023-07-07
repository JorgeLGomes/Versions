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
export Dir_home=/home/sondamod/CFS
export Dir_wrk=/rede/etap01/sonda/CFS
export Dir_wrk_out=${Dir_wrk}/data/CFS/output/${Run_Date}
export Dir_scr=${Dir_home}/scripts
export Dir_util=${Dir_home}/util
export data_grib=${Dir_wrk}/data/CFS
export Fct=264
export FInitBC=CFS
export InitBC=6

ifct=000
ihour=`echo ${Run_Date} |cut -c 9-10`


export GRBDIR=${data_grib}/${Run_Date}
export WOUTDIR=${Dir_wrk_out}
mkdir -p ${Dir_wrk_out}
me=$(whoami)
echo "I am $me."


dhri=${ifct}

typeset -Z4 dhr
typeset -Z3 dhri Fct
while [ $dhri -le ${Fct} ];do
echo "passei"
xargs -P 8 |CFS_deco_parallel.ksh $1 $dhri
dhri=`expr $dhri + ${InitBC} `

done # while [ $dhr -le $dend ]


exit
