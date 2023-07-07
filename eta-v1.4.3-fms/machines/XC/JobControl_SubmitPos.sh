 NJobs=`qstat |grep ${Exp}${Run_Date_Hex} | grep ${USER} | wc -l`
 if ((${NJobs}<7)) ; then
