#!/bin/bash -x
#
# Verifica se foi passado algum parametro
if (($#<2)) ; then
echo "Use: gfs2_deco.sh yyyymmddhh fct"
fi
Run_Date=${1}
ifct=${2}

# VARIAVEIS
export Dir_scr=DIRROOT
export Dir_home=`dirname ${Dir_scr}` 
export Dir_datain=`dirname ${Dir_home}` 
export ModelDrive=`basename  ${Dir_scr}`
export Dir_wrk=${Dir_datain}/atmos/${ModelDrive}/${Run_Date}
export Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
export Dir_util=${Dir_datain}/util
export Dir_dprep_exe=${Dir_datain}/dprep/exe
export InitBC=6
export IntSST=24
export FInitBC=${ModelDrive}



cp ${Dir_scr}/InputModelInf.txt_${FInitBC} ${Dir_ETAwrk}/InputModelInf.txt
dprep_exe=`head -1 ${Dir_ETAwrk}/InputModelInf.txt`

me=$(whoami)
echo "I am $me."

dhri=${ifct}
dhriF=`printf "%03d" "${dhri}"`

ihour=`echo ${Run_Date} |cut -c 9-10`

dhr=`printf "%06d" "${dhri}"`

ifile='gfs.t'${ihour}'z.pgrb2.0p25.f'${dhriF}
echo '==================================================='
echo 'file in process '${Dir_wrk}'/'${ifile}
echo '==================================================='

line=$(${Dir_util}/wgrib2 -v ${Dir_wrk}'/'${ifile} | grep "HGT" | grep "10 mb")
adate=`echo ${line}|cut -d"=" -f2|cut -c1-10`
echo $adate > ${Dir_ETAwrk}/adate.txt

if [ $adate -eq $Run_Date ] ; then
  echo "The gfs date is $adate"
else
  echo "Data for $Run_Date isn't available"
  exit 99 
fi

rm -f  ${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr
touch  ${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr
echo 'OPEN ' ${Dir_wrk}/${FInitBC}_${adate}.$dhr '  FILE'
if [[ ! -s ${Dir_ETAwrk}/gfs2_field_rec.txt ]] ; then
for vars in HGT  RH  SPFH TMP UGRD VGRD; do
for levs in 100000 97500 95000 92500 90000 85000 80000 75000 70000 65000 \
             60000 55000 50000 45000 40000 35000 30000 25000 20000 15000 \
             10000  7000  5000  3000  2000  1000;do


line=$(${Dir_util}/wgrib2 -ctl_inv ${Dir_wrk}'/'${ifile} | grep "$vars 100,$levs ")

if [ -z "$line" ]; then
echo 'MISSING DATA :: VAR = '$vars' LEVEL = '$levs
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec ${Dir_wrk}'/'${ifile}  -append -bin ${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr >> ${Dir_ETAwrk}/log.${dhriF}
echo  "${nrec}" >> ${Dir_ETAwrk}/gfs2_field_rec.txt
fi   # if [ -z "$line" ]

done # for levs 
done # for vars 

# SURFACE VARIABLES

for vars in "LAND" "PRMSL" "PRES 1," "HGT 1,"\
            "TSOIL 106,0,0.1"\
            "SOILW 106,0,0.1"\
            "TSOIL 106,0.1,0.4"\
            "SOILW 106,0.1.0.4"\
            "TSOIL 106,0.4,1"\
            "SOILW 106,0.4,1"\
            "TSOIL 106,1,2"\
            "SOILW 106,1,2"\
	    ; do

line=$(${Dir_util}/wgrib2 -ctl_inv ${Dir_wrk}'/'${ifile} | grep "$vars")

if [ -z "$line" ]; then
echo 'NO DATA FOR VAR = '$vars' LEVEL = '$levs'mb'
else

tptp=`expr index "$line" :`
tptp=`expr $tptp - 1`
nrec=${line:0:$tptp}

${Dir_util}/wgrib2 -d $nrec ${Dir_wrk}'/'${ifile}  -append -bin  ${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr >> ${Dir_ETAwrk}/log.${dhriF}
echo  "${nrec}" >> ${Dir_ETAwrk}/gfs2_field_rec.txt
fi   # if [ -z "$line" ]
done # for vars


echo 'DONE WITH ' ${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr '  FILE'
else

cat ${Dir_ETAwrk}/gfs2_field_rec.txt | while read nrec
do
   ${Dir_util}/wgrib2 -d $nrec ${Dir_wrk}'/'${ifile}  -append -bin ${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr >> ${Dir_ETAwrk}/log.${dhriF}
   echo "nrec: $nrec"
done
echo 'DONE WITH ' ${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr '  FILE'
ls -l ${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr
fi
ExecCmd1P ${Dir_dprep_exe}/${dprep_exe} << endin
${Dir_ETAwrk}/${FInitBC}_${adate}.$dhr
${Dir_ETAwrk}
endin
exit
