#!/bin/bash -x
#
if (($#<2)) ; then
   echo 'Parameters: yyyymmddhh (init) yyyymmddhh (end)'   
   exit
fi
Run_Date=${1}
datef=${2}
# VARIAVEIS
Dir_scr=DIRROOT
Dir_home=`dirname ${Dir_scr}`
Dir_datain=`dirname ${Dir_home}`
ModelDrive=`basename  ${Dir_scr}`
Dir_wrk=${Dir_datain}/atmos/${ModelDrive}/${Run_Date}
Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
Dir_SST=${Dir_datain}/sst/${ModelDrive}
Dir_util=${Dir_datain}/util
InitBC=6
FInitBC=${ModelDrive}
mkdir -p ${Dir_ETAwrk}

let FctEnd=(`date +%s -d ${datef:0:8}`-`date +%s -d ${Run_Date:0:8}`)/3600
echo ${FctEnd}
rm -f  ${Dir_ETAwrk}/Submit_process.list
LastArq="False"
tval=0
tvalF=`printf "%06d" "${tval}"`
FctEndF=`printf "%06d" "${FctEnd}"`
while [ "${LastArq}" != "True" ] ; do
cat <<EOF> ${Dir_ETAwrk}/Submit_process${tvalF}
#!/bin/bash
cd ${Dir_ETAwrk}
${Dir_scr}/process.sh ${Run_Date} ${tval}
EOF
chmod 755 ${Dir_ETAwrk}/Submit_process${tvalF}
echo "${Dir_ETAwrk}/Submit_process${tvalF}" >> ${Dir_ETAwrk}/Submit_process.list
if [ "${tvalF}" == "${FctEndF}" ] ; then
  GlobalOK="True"
  break
fi
  let tval=${tval}+${InitBC}
  tvalF=`printf "%06d" "${tval}"`
done
cat ${Dir_ETAwrk}/Submit_process.list  | xargs -n 1 -P 8 /bin/bash
#limpeza
rm -f ${Dir_ETAwrk}/Submit_process??????
rm -f ${Dir_ETAwrk}/${FInitBC}_${Run_Date}.??????
rm -f ${Dir_ETAwrk}/Submit_process.list
exit
