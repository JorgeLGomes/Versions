#!/bin/ksh -x
#
# Verifica se foi passado algum parametro
if (($#<4)) ; then
export Run_Date=`date "+%Y%m%d"`$1
else
export Run_Date=${4}${1}
fi

export CFS_Member=${2}
export Fcti=${3}

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

cp ${Dir_scr}/InputModelInf.txt_${ModelDrive} ${Dir_ETAwrk}/InputModelInf.txt

me=$(whoami)
echo "I am $me."

echo 'remove old files!!!'
rm -fv wbin.*
rm -fv wbin.*


typeset -Z6 Fcti 

datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${Fcti}hr 'yyyymmddhh'`
ifilepgbf=${Dir_wrk}'/pgbf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2'
ifileflxf=${Dir_wrk}'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2'
ofile=${Dir_ETAwrk}'/'${FInitBC}'_'${Run_Date}'.'${Fcti}
echo '==================================================='
echo 'file in process '${ifilepgbf}
echo '==================================================='

if [ ${Fcti} -eq 000000 ];then
  if [[ ! -s ${ifileflxf} ]] ; then
    datefct2=`${Dir_util}/caldate.3.0 ${Run_Date} + 6hr 'yyyymmddhh'` 
    ifileflxf=${Dir_wrk}'/flxf'${datefct2}'.'${CFS_Member}'.'${Run_Date}'.grb2'
  fi
line=$(${Dir_util}/wgrib2 -v ${ifilepgbf} | grep "HGT" | grep "10 mb")
adate=`echo ${line}|cut -d"=" -f2|cut -c1-10`
echo $adate > ${Dir_ETAwrk}/adate.txt
echo $adate $Run_Date

if [ $adate -eq $Run_Date ] ; then
  echo "The cfs date is $adate"
else
  echo "Data for $Run_Date isn't available"
  exit 99 
fi
fi

rm -f  ${ofile}
touch  ${ofile}
echo 'OPEN ' ${ofile} '  FILE'

if [[ ! -s ${Dir_ETAwrk}/cfs_field_rec.txt ]] ; then
for vars in HGT RH SPFH TMP UGRD VGRD; do
for levs in 100000 97500 95000 92500 90000 85000 80000 75000 70000 65000 \
             60000 55000 50000 45000 40000 35000 30000 25000 20000 15000 \
             10000  7000  5000  3000  2000  1000;do


line=$(${Dir_util}/wgrib2 -ctl_inv ${ifilepgbf} | grep "$vars 100,$levs ")

if [ -z "$line" ]; then
echo 'MISSING DATA :: VAR = '$vars' LEVEL = '$levs
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}


${Dir_util}/wgrib2 -d $nrec ${ifilepgbf} -ncpu 4 -append -bin ${ofile} >> ${Dir_ETAwrk}/log.${Fcti}
echo  "${nrec}" >> ${Dir_ETAwrk}/cfs_field_rec.txt
fi   # if [ -z "$line" ]

done # for levs 
done # for vars 
else
cat ${Dir_ETAwrk}/cfs_field_rec.txt | while read nrec
do
   ${Dir_util}/wgrib2 -d $nrec ${ifilepgbf} \
   -ncpu 4 -append -bin ${ofile} >> ${Dir_ETAwrk}/log.${Fcti}
   echo "nrec: $nrec"
done

fi 
# SURFACE VARIABLES

# LAND

line=$(${Dir_util}/wgrib2 -ctl_inv ${ifileflxf} | grep "LAND")

if [ -z "$line" ]; then
echo 'NO DATA FOR LAND'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec ${ifileflxf} -ncpu 4 -new_grid latlon 0:360:1 -90:181:1 ${Dir_ETAwrk}/append.grb2.${Fcti} >> ${Dir_ETAwrk}/log.${Fcti}
${Dir_util}/wgrib2 -d 1 ${Dir_ETAwrk}/append.grb2.${Fcti}  -append -bin ${ofile} >> ${Dir_ETAwrk}/log.${Fcti}

fi

# PRMSL

line=$(${Dir_util}/wgrib2 -ctl_inv ${ifilepgbf} | grep "PRMSL")

if [ -z "$line" ]; then

echo 'NO DATA FOR PRMSL'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec ${ifilepgbf} -ncpu 4 -append -bin  ${ofile} >> ${Dir_ETAwrk}/log.${Fcti}

fi

# PRESsfc

line=$(${Dir_util}/wgrib2 -ctl_inv ${ifileflxf} | grep "PRES 1")

if [ -z "$line" ]; then

echo 'NO DATA FOR PRES 1'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}


${Dir_util}/wgrib2 -d $nrec ${ifileflxf} -ncpu 4 -new_grid latlon 0:360:1 -90:181:1  ${Dir_ETAwrk}/append.grb2.${Fcti} >> ${Dir_ETAwrk}/log.${Fcti}
${Dir_util}/wgrib2 -d 1 ${Dir_ETAwrk}/append.grb2.${Fcti}  -ncpu 4 -append -bin ${ofile} >> ${Dir_ETAwrk}/log.${Fcti}
fi

# HGTsfc

line=$(${Dir_util}/wgrib2 -ctl_inv ${ifileflxf} | grep "HGT 1")

if [ -z "$line" ]; then

echo 'NO DATA FOR HGTsfc'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec ${ifileflxf} -ncpu 4 -new_grid latlon 0:360:1 -90:181:1  ${Dir_ETAwrk}/append.grb2.${Fcti} >> ${Dir_ETAwrk}/log.${Fcti}
${Dir_util}/wgrib2 -d 1 ${Dir_ETAwrk}/append.grb2.${Fcti}  -ncpu 4 -append -bin ${ofile} >> ${Dir_ETAwrk}/log.${Fcti}

fi

# SOIL VARIABLES

if [[ ! -s ${Dir_ETAwrk}/cfs_soil_rec.txt ]] ; then
for vars in "TMP 106,0,0.1"\
            "SOILW 106,0,0.1"\
            "TMP 106,0.1,0.4"\
            "SOILW 106,0.1.0.4"\
            "TMP 106,0.4,1"\
            "SOILW 106,0.4,1"\
            "TMP 106,1,2"\
            "SOILW 106,1,2"\
            ; do

line=$(${Dir_util}/wgrib2 -ctl_inv ${ifileflxf} | grep "$vars")

if [ -z "$line" ]; then
echo 'NO DATA FOR VAR = '$vars' LEVEL = '$levs'mb'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec ${ifileflxf} -ncpu 4 -new_grid latlon 0:360:1 -90:181:1  ${Dir_ETAwrk}/append.grb2.${Fcti} >> ${Dir_ETAwrk}/log.${Fcti}
${Dir_util}/wgrib2 -d 1 ${Dir_ETAwrk}/append.grb2.${Fcti}  -ncpu 4 -append -bin ${ofile} >> ${Dir_ETAwrk}/log.${Fcti}
echo  "${nrec}" >> ${Dir_ETAwrk}/cfs_soil_rec.txt
rm -f ${Dir_ETAwrk}/append.grb2.${Fcti}
fi   # if [ -z "$line" ]
done # for vars

else
cat ${Dir_ETAwrk}/cfs_soil_rec.txt | while read nrec
do
   ${Dir_util}/wgrib2 -d $nrec ${ifileflxf} -ncpu 4 -new_grid latlon 0:360:1 -90:181:1  ${Dir_ETAwrk}/append.grb2.${Fcti}  >> ${Dir_ETAwrk}/log.${Fcti}
   ${Dir_util}/wgrib2 -d 1 ${Dir_ETAwrk}/append.grb2.${Fcti} -ncpu 4 -append -bin ${ofile} >> ${Dir_ETAwrk}/log.${Fcti}
   echo "nrec: $nrec"
   rm -f ${Dir_ETAwrk}/append.grb2.${Fcti}
done

fi 

echo 'DONE WITH ' ${ofile} '  FILE'

ls -l ${ofile}
ExecCmd1P ${Dir_dprep_exe}/dg${ModelDrive}.exe << endin
${ofile}
${Dir_ETAwrk}
endin

if [ -s ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}.ETAwrk ]  ; then
  rm -f ${ofile}
  exit 0
else
ExecCmd1P ${Dir_dprep_exe}/dg${ModelDrive}.exe << endin
${ofile}
${Dir_ETAwrk}
endin
echo "${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}.ETAwrk" >> ${Dir_ETAwrk}/${ModelDrive}.ETAwrk.list
  if [ ! -s ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}.ETAwrk ]  ; then
    echo "Verificar ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}" >> ${Dir_ETAwrk}/${ModelDrive}.ETAwrk.list
    exit 1
  fi
fi

exit
