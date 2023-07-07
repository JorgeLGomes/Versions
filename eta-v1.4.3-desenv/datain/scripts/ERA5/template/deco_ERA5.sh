#!/bin/ksh -x
#
# Verifica se foi passado algum parametro

export Run_Date=${1}


# VARIAVEIS
Dir_scr=DIRROOT
Dir_home=`dirname ${Dir_scr}`
Dir_datain=`dirname ${Dir_home}`
ModelDrive=`basename ${Dir_scr}`
Dir_wrk=${Dir_datain}/atmos/${ModelDrive}/${Run_Date}
Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
Dir_ETAwrk_soil=${Dir_datain}/soil/ETAwrk/${ModelDrive}/${Run_Date}
Dir_util=${Dir_datain}/util
Dir_dprep_exe=${Dir_datain}/dprep/exe
InitBC=6
FInitBC=${ModelDrive}

if [ ! -d ${Dir_ETAwrk} ] ; then
  mkdir -p  ${Dir_ETAwrk}
fi
if [ ! -d ${Dir_ETAwrk_soil} ] ; then
  mkdir -p  ${Dir_ETAwrk_soil}
fi
cd ${Dir_ETAwrk}
cp ${Dir_scr}/InputModelInf.txt_${ModelDrive} ${Dir_ETAwrk}/InputModelInf.txt

me=$(whoami)
echo "I am $me."

echo 'remove old files!!!'
rm -fv wbin.*
rm -fv wbin.*


typeset -Z6 Fcti 

Fcti=000000

echo '==================================================='
echo 'file in process '${Dir_wrk}'/ERA5_Pressure_'${Run_Date}'.grib'
echo '==================================================='

line=$(${Dir_util}/wgrib -v ${Dir_wrk}'/ERA5_Pressure_'${Run_Date}'.grib' | grep "Geopotential" | grep "20 mb")
adate=`echo ${line}|cut -d"=" -f2|cut -c1-10`
echo $adate > ${Dir_ETAwrk}/adate.txt
echo $adate $Run_Date

if [ $adate -eq $Run_Date ] ; then
  echo "The ERA5 date is $adate"
else
  echo "Data for $Run_Date isn't available"
  exit 99 
fi

rm -f  ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}
touch  ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}
echo 'OPEN ' ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti} '  FILE'

nrec=1
while ((${nrec}<=222)); do
   ${Dir_util}/wgrib -d $nrec ${Dir_wrk}'/ERA5_Pressure_'${Run_Date}'.grib' \
                     -nh -bin  -append -o ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti} >> ${Dir_ETAwrk}/log.${Fcti}
let nrec=${nrec}+1
done
nrec=1
while ((${nrec}<=12)); do
   ${Dir_util}/wgrib -d $nrec ${Dir_wrk}'/ERA5_Surface_'${Run_Date}'.grib' \
                     -nh -bin  -append -o ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti} >> ${Dir_ETAwrk}/log.${Fcti}
let nrec=${nrec}+1
done

cd ${Dir_ETAwrk_soil}
cp ${Dir_scr}/InputModelInf.txt_${ModelDrive} ${Dir_ETAwrk_soil}/InputModelInf.txt
nrec=5
while ((${nrec}<=12)); do
   ${Dir_util}/wgrib -d $nrec ${Dir_wrk}'/ERA5_Surface_'${Run_Date}'.grib' \
                     -nh -bin  -append -o ${Dir_ETAwrk_soil}/${FInitBC}_SoilT_${Run_Date}.${Fcti} >> ${Dir_ETAwrk_soil}/log.${Fcti}
let nrec=${nrec}+2
done
nrec=6
while ((${nrec}<=12)); do
   ${Dir_util}/wgrib -d $nrec ${Dir_wrk}'/ERA5_Surface_'${Run_Date}'.grib' \
                     -nh -bin  -append -o ${Dir_ETAwrk_soil}/${FInitBC}_SoilW_${Run_Date}.${Fcti} >> ${Dir_ETAwrk_soil}/log.${Fcti}
let nrec=${nrec}+2
done

echo 'DONE WITH ' ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti} '  FILE'

cd ${Dir_ETAwrk}
cp ${Dir_scr}/InputModelInf.txt_${ModelDrive} ${Dir_ETAwrk}/InputModelInf.txt
ls -l ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}
ExecCmd1P ${Dir_dprep_exe}/dg${ModelDrive}.exe << endin
${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}
${Dir_ETAwrk}
endin
rm -f ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.${Fcti}
rm -f ${Dir_wrk}/ERA5_Pressure_${Run_Date}.grib
rm -f ${Dir_wrk}/ERA5_Surface_${Run_Date}.grib
rm -f ${Dir_wrk}/Submit_deco${Run_Date}
rm -f ${Dir_ETAwrk}/log.${Fcti}

exit
