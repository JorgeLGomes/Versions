#!/bin/ksh -x
#
# Verifica se foi passado algum parametro

export Run_Date=${1}


# VARIAVEIS
Dir_scr=DIRROOT
Dir_home=`dirname ${Dir_scr}`
Dir_datain=`dirname ${Dir_home}`
ModelDrive=`basename ${Dir_scr}`
Dir_wrk=${Dir_datain}/sst/${ModelDrive}/${Run_Date}
Dir_ETAwrk=${Dir_datain}/sst/ETAwrk/${ModelDrive}/${Run_Date}
Dir_util=${Dir_datain}/util
InitBC=6
FInitBC=${ModelDrive}

if [ ! -d ${Dir_ETAwrk} ] ; then
  mkdir -p  ${Dir_ETAwrk}
fi
cd ${Dir_ETAwrk}

me=$(whoami)
echo "I am $me."

typeset -Z6 Fcti 
typeset -Z4 TNTimes

TNTimes=1

Fcti=000000
ifile=${Dir_ETAwrk}/ERA5_SST_${Run_Date}.grib
ofile=${Dir_ETAwrk}/${ModelDrive}_sst_${Run_Date}.bin
echo '==================================================='
echo 'file in process '${ifile}
echo '==================================================='
rm -f  ${ofile}
${Dir_util}/wgrib -d 1 ${ifile} \
                  -nh -bin  -append -o ${ofile} >> ${Dir_ETAwrk}/log.${Run_Date}

echo 'DONE WITH ' ${ofile}'  FILE'
echo ${TNTimes} >${Dir_ETAwrk}/NTimes.txt
ls -l ${ofile}

rm -f ${ifile}

exit
