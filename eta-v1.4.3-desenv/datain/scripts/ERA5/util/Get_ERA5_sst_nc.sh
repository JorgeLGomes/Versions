#!/bin/bash 
#
# Get_ERA5_nc+process.sh yyyymmddhh nhour
# yyyymmddhh : initial date to start to get data
#       nhour: number of hours of data download 
#              frequency 6/6h

# Check input parameter 
if (($#<2)) ; then
  echo "      Get_ERA5_sst_nc.sh yyyymmddhh nhour"
  echo "      yyyymmddhh : initial date to start to get data"
  echo "      nhour: number of hours of data download"
  echo "      frequency 24/24h"
  exit 1
fi

export Run_Date=${1}
export nhours=${2}
yyyyi=${Run_Date:0:4} 
mmi=${Run_Date:4:2}
ddi=${Run_Date:6:2} 
hhi=${Run_Date:10:2}

# VARIAVEIS
Dir_scr=/home/jorge/EtaLab/eta/datain/scripts/ERA5
Dir_home=`dirname ${Dir_scr}`
Dir_datain=`dirname ${Dir_home}`
ModelDrive=`basename ${Dir_scr}`
Dir_wrk=${Dir_datain}/atmos/${ModelDrive}/${Run_Date}
Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
Dir_ETAwrk_SST=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
Dir_util=${Dir_datain}/util
InitBC=24
FInitBC=${ModelDrive}

LastArq="False"
tval=000000
tvalF=`printf "%06d" "${tval}"`
nhours=`printf "%06d" "${nhours}"`


rm -f ${Dir_scr}/Submit_deco.list

while [ "${LastArq}" != "True" ] ; do
Dir_wrk=${Dir_datain}/atmos/${ModelDrive}/${Run_Date}
Dir_ETAwrk=${Dir_datain}/atmos/ETAwrk/${ModelDrive}/${Run_Date}
Dir_ETAwrk_SST=${Dir_datain}/sst/ETAwrk/${ModelDrive}/${Run_Date}
if [ ! -d ${Dir_wrk} ] ; then
  mkdir -p  ${Dir_wrk}
fi
if [ ! -d ${Dir_ETAwrk} ] ; then
  mkdir -p  ${Dir_ETAwrk}
fi
if [ ! -d ${Dir_ETAwrk_SST} ] ; then
  mkdir -p  ${Dir_ETAwrk_SST}
fi
cat <<EOF> ${Dir_wrk}/Submit_deco${Run_Date}
#!/bin/bash -x

cd ${Dir_wrk}
python3 ${Dir_scr}/get_ERA5_sst_nc.py ${Run_Date:0:4} ${Run_Date:4:2} ${Run_Date:6:2} ${Run_Date:8:2} 
EOF
  chmod 755 ${Dir_wrk}/Submit_deco${Run_Date}
  echo "${Dir_wrk}/Submit_deco${Run_Date}" >> ${Dir_scr}/Submit_deco.list
  if [ "${tvalF}" == "${nhours}" ] ; then
    GlobalOK="True"
    break
  fi
  let tval=${tval}+${InitBC}
  tvalF=`printf "%06d" "${tval}"`
  Run_Date=`date -d "${yyyyi}-${mmi}-${ddi} ${hhi} +${tval} hour" "+%Y%m%d%H"`
done

cat ${Dir_scr}/Submit_deco.list  | xargs -n 1 -P 1 /bin/bash
#limpeza
rm -f ${Dir_ETAwrk}/Submit_deco.list
exit
