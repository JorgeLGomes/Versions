#!/bin/bash -x
#
# Verifica se foi passado algum parametro
if (($#<2)) ; then
echo "Use: "
fi
Run_Date=${1}
ifct=${2}

# VARIAVEIS
export Dir_scr=DIRROOT
export Dir_home=`dirname ${Dir_scr}` 
export Dir_datain=`dirname ${Dir_home}` 
export ModelDrive=`basename ${Dir_home}`
export Dir_wrk=${Dir_datain}/atmos/${ModelDrive}/${Run_Date}
export Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
export Dir_util=${Dir_datain}/util
export Dir_dprep_exe=${Dir_datain}/dprep/exe
export InitBC=6
export IntSST=24
export FInitBC=gposeta

FInitBC_Up=`echo ${FInitBC}|awk '{print toupper($1)}'`
cp ${Dir_scr}/InputModelInf.txt_${FInitBC} ${Dir_ETAwrk}/InputModelInf.txt

me=$(whoami)
echo "I am $me."

cd ${Dir_wrk}

dhr=`printf "%06d" "${ifct}"`
if [ "${dhr}" == "000000" ]  ; then
   datef=${Run_Date}
   ffct=icn
else
   datei=${Run_Date}
   datef=`${Dir_util}/caldate.3.0 ${Run_Date}   + ${InitBC}hr 'yyyymmddhh'`
   ffct=fct
fi
ifile=${Dir_wrk}/${FInitBC_Up}${Run_Date}${datef}E.${ffct}.TQ0062L042
ofile=${Dir_ETAwrk}/${ModelDrive}${Run_Date}.${dhr}
echo '==================================================='
echo 'file in process '${ifile}
echo '==================================================='

echo ${Run_Date} > ${Dir_ETAwrk}/adate.txt
ln -s ${ifile} ${ofile}
${Dir_dprep_exe}/dg${FInitBC}.exe << endin
${ofile}
${Dir_ETAwrk}
endin
exit
