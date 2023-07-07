#!/bin/ksh -x
#
# Verifica se foi passado algum parametro
if (($#<5)) ; then
export Run_Date=`date "+"%Y%m%d"`$1
else
export Run_Date=${5}${1}
fi

# VARIAVEIS
export FInitBC=CFS
export CFS_Member=${2}
export Dir_home=/dados/grpeta/dsk003/scripts/CFS
export Dir_wrk=/dados/grpeta/dsk003/dados/CFS
export Dir_SST=/dados/grpeta/dsk003/dados/SST
export Dir_scr=${Dir_home}/scripts
export Dir_wrk_out=${Dir_wrk}/output/${Run_Date}_M${CFS_Member}
export Dir_util=${Dir_home}/util
export data_grib=${Dir_wrk}
export ifct=${3}
export Fct=${4}
export InitBC=6

#ifct=0000
ihour=`echo ${Run_Date} |cut -c 9-10`


export GRBDIR=${data_grib}/${Run_Date}_M${CFS_Member}
export WOUTDIR=${Dir_wrk_out}
mkdir -p ${Dir_wrk_out}
me=$(whoami)
echo "I am $me."

echo 'remove old files!!!'
rm -fv $WOUTDIR/wbin.*
rm -fv wbin.*

dhri=${ifct}

typeset -Z4 dhr Fct
typeset -Z4 dhri 
while [ $dhri -le ${Fct} ];do
echo "passei"

#if [ $dhri -ne 0 ];then
#if [ $dhri -le 9 ];then
#dhri=0$dhri
#fi
#fi

dhr=${dhri}
datefct=`${Dir_util}/caldate.3.0 ${Run_Date} + ${dhri}hr 'yyyymmddhh'`

echo '==================================================='
echo 'file in process '$GRBDIR'/pgbf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2'
echo '==================================================='

#if [ $dhri -eq 0 ];then

line=$(${Dir_util}/wgrib2 -v $GRBDIR'/pgbf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "HGT" | grep "10 mb")
adate=`echo ${line}|cut -d"=" -f2|cut -c1-10`
echo $adate > ${WOUTDIR}/adate.txt
echo $adate $Run_Date

if [ $adate -eq $Run_Date ] ; then
  echo "The cfs date is $adate"
else
  echo "Data for $Run_Date isn't available"
  exit 99 
fi
#fi

rm -f  ${WOUTDIR}/${FInitBC}_${adate}E.$dhr
touch  ${WOUTDIR}/${FInitBC}_${adate}E.$dhr
echo 'OPEN ' ${WOUTDIR}/${FInitBC}_${adate}E.$dhr '  FILE'


for vars in HGT UGRD VGRD TMP RH; do
for levs in 100000 97500 95000 92500 90000 85000 80000 75000 70000 65000 \
             60000 55000 50000 45000 40000 35000 30000 25000 20000 15000 \
             10000  7000  5000  3000  2000  1000;do


line=$(${Dir_util}/wgrib2 -ctl_inv $GRBDIR'/pgbf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "$vars 100,$levs ")

if [ -z "$line" ]; then
echo 'MISSING DATA :: VAR = '$vars' LEVEL = '$levs
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}


${Dir_util}/wgrib2 -d $nrec $GRBDIR'/pgbf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' -append -bin ${WOUTDIR}/${FInitBC}_${adate}E.$dhr >> ${WOUTDIR}/log

fi   # if [ -z "$line" ]

done # for levs 
done # for vars 

# SURFACE VARIABLES

# LAND

line=$(${Dir_util}/wgrib2 -ctl_inv $GRBDIR'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "LAND")

if [ -z "$line" ]; then
echo 'NO DATA FOR LAND'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec $GRBDIR'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' -new_grid latlon 0:360:1 -90:181:1 ${WOUTDIR}/append.grb2 >> ${WOUTDIR}/log
${Dir_util}/wgrib2 -d 1 ${WOUTDIR}/append.grb2  -append -bin ${WOUTDIR}/${FInitBC}_${adate}E.$dhr >> ${WOUTDIR}/log

fi

# PRMSL

line=$(${Dir_util}/wgrib2 -ctl_inv $GRBDIR'/pgbf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "PRMSL")

if [ -z "$line" ]; then

echo 'NO DATA FOR PRMSL'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec $GRBDIR'/pgbf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' -append -bin  ${WOUTDIR}/${FInitBC}_${adate}E.$dhr >> ${WOUTDIR}/log

fi

# PRESsfc

line=$(${Dir_util}/wgrib2 -ctl_inv $GRBDIR'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "PRES 1")

if [ -z "$line" ]; then

echo 'NO DATA FOR PRES 1'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}


${Dir_util}/wgrib2 -d $nrec $GRBDIR'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' -new_grid latlon 0:360:1 -90:181:1  ${WOUTDIR}/append.grb2 >> ${WOUTDIR}/log
${Dir_util}/wgrib2 -d 1 ${WOUTDIR}/append.grb2  -append -bin ${WOUTDIR}/${FInitBC}_${adate}E.$dhr >> ${WOUTDIR}/log
fi

# HGTsfc

line=$(${Dir_util}/wgrib2 -ctl_inv $GRBDIR'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "HGT 1")

if [ -z "$line" ]; then

echo 'NO DATA FOR HGTsfc'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec $GRBDIR'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' -new_grid latlon 0:360:1 -90:181:1  ${WOUTDIR}/append.grb2 >> ${WOUTDIR}/log
${Dir_util}/wgrib2 -d 1 ${WOUTDIR}/append.grb2  -append -bin ${WOUTDIR}/${FInitBC}_${adate}E.$dhr >> ${WOUTDIR}/log

fi


# SOIL VARIABLES

for vars in "TMP 106,0,0.1"\
            "SOILW 106,0,0.1"\
            "TMP 106,0.1,0.4"\
            "SOILW 106,0.1.0.4"\
            "TMP 106,0.4,1"\
            "SOILW 106,0.4,1"\
            "TMP 106,1,2"\
            "SOILW 106,1,2"\
            ; do

line=$(${Dir_util}/wgrib2 -ctl_inv $GRBDIR'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' | grep "$vars")

if [ -z "$line" ]; then
echo 'NO DATA FOR VAR = '$vars' LEVEL = '$levs'mb'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec $GRBDIR'/flxf'${datefct}'.'${CFS_Member}'.'${Run_Date}'.grb2' -new_grid latlon 0:360:1 -90:181:1  ${WOUTDIR}/append.grb2 >> ${WOUTDIR}/log
${Dir_util}/wgrib2 -d 1 ${WOUTDIR}/append.grb2  -append -bin ${WOUTDIR}/${FInitBC}_${adate}E.$dhr >> ${WOUTDIR}/log

fi   # if [ -z "$line" ]
done # for vars




echo 'DONE WITH ' ${WOUTDIR}/${FInitBC}_${adate}E.$dhr '  FILE'



ls -l ${WOUTDIR}/${FInitBC}_${adate}E.$dhr

#echo '=========================================='
dhri=`expr $dhri + ${InitBC} `

done # while [ $dhr -le $dend ]


exit
