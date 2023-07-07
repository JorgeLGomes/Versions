#!/bin/bash -x
#
#
# Get+process.sh hh member fcst_ini fcst_end yyyymmdd
# 
# Check input parameter
if (($#<5)) ; then
  echo " Get+process.sh hh member fcst_ini fcst_end yyyymmdd"
  echo "       hh: initial condition hour"
  echo "   member: CFS member (this script gets only control"
  echo "           member 01"
  echo " fcst_ini: initial forecast => 00"
  echo " fcst_end: number of hours of data download"
  echo " yyyymmdd: initial condition date of data"
  exit 1
fi

hostname
export hh=${1}
if (($#<5)) ; then
export Run_Date=`date "+%Y%m%d"`${hh}
else
export Run_Date=${5}${hh}
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
export InitBC=6
export IntSST=24


mkdir -p ${Dir_wrk}
mkdir -p ${Dir_ETAwrk}
mkdir -p ${Dir_SST_ETAwrk}

${Dir_scr}/get_CFS_wgetParallel.ksh ${hh} ${CFS_Member} ${Fcti} ${Fct} ${Run_Date:0:8}

${Dir_scr}/sstcfs_grb2_deco.sh ${hh}  ${CFS_Member} 24 ${Fct} ${Run_Date:0:8}
rm -f  ${Dir_ETAwrk}/Submit_deco.list

# PbsSlurm

LastArq="False"
tval=`printf "%06d" "${Fcti}"`
Fct=`printf "%06d" "${Fct}"`
yymmdd=`echo ${Run_Date} |cut -c3-8`
while [ "${LastArq}" != "True" ] ; do
cat <<EOF> ${Dir_ETAwrk}/Submit_deco${tval}
#!/bin/bash
cd ${Dir_ETAwrk}
${Dir_scr}/CFS_deco.ksh ${hh}  ${CFS_Member} ${tval} ${Run_Date:0:8}
EOF
  chmod 755 ${Dir_ETAwrk}/Submit_deco${tval}
  echo "${Dir_ETAwrk}/Submit_deco${tval}" >> ${Dir_ETAwrk}/Submit_deco.list
  if [ 10#${tval} == 10#${Fct} ] ; then
    GlobalOK="True"
    break
  fi
  let tval=10#${tval}+${InitBC}
  tval=`printf "%06d" "${tval}"`
done

# QueueCmd

#limpeza
rm -f ${Dir_ETAwrk}/Submit_deco??????
rm -f  ${Dir_ETAwrk}/log.???
rm -f  ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.??????
rm -f ${Dir_ETAwrk}/Submit_deco.list
rm -f ${Dir_ETAwrk}/cfs_field_rec.txt
rm -f ${Dir_ETAwrk}/cfs_soil_rec.txt
exit
