#!/bin/bash -x

LonW=$1
LonE=$2
LatN=$3
LatS=$4
res=$5
#INCLUI os diretorios de scripts e de saida
export LC_NUMERIC="en_US.UTF-8"
Eta_ucl=../ucl

west1=`echo ${LonW} | awk '{printf "%11.6f",$1}'`
#echo $west1
west2=`echo ${LonW} |awk 'invertsgn = $1 * (-1) {printf "%11.6f",invertsgn}'`
#echo $west2
south=`echo ${LatS} | awk '{printf "%11.6f",$1}'`
#echo $south
res=`echo $res| awk '{printf "%11.6f",$1}'`
#echo $res
imout=`expr " (((${LonE} - ${LonW})  / $res ) + 1 )" |bc -l| awk '{printf "%5.0f",$1}'`
imout1=`printf "%05d" "${imout}"`
jmout=`expr " (((${LatN} - ${LatS})  / $res ) + 1 )" |bc -l| awk '{printf "%5.0f",$1}'`
jmout1=`printf "%05d" "${jmout}"`

cat ${Eta_ucl}/CTLTEMPLATE_2D_temp    | sed "s:ptx:$imout1:g"       \
                                      | sed "s:pty:$jmout1:g"       \
                                      | sed "s:lonoeste:${west1}:g" \
                                      | sed "s:latsul:${south}:g"   \
  > ${Eta_ucl}/CTLTEMPLATE_2D
#
cat ${Eta_ucl}/CTLTEMPLATE_3D_temp    | sed "s:ptx:$imout1:g"       \
                                      | sed "s:pty:$jmout1:g"       \
                                      | sed "s:lonoeste:${west1}:g" \
                                      | sed "s:latsul:${south}:g"   \
  > ${Eta_ucl}/CTLTEMPLATE_3D
#
cat ${Eta_ucl}/CTLTEMPLATE_2D3D_temp  | sed "s:ptx:$imout1:g"       \
                                      |  sed "s:pty:$jmout1:g"      \
                                      | sed "s:lonoeste:${west1}:g" \
                                      | sed "s:latsul:${south}:g"   \
  > ${Eta_ucl}/CTLTEMPLATE_2D3D
#
cat ${Eta_ucl}/CTLTEMPLATE_FIXED_temp | sed "s:ptx:$imout1:g"       \
                                      | sed "s:pty:$jmout1:g"       \
                                      | sed "s:lonoeste:${west1}:g" \
                                      | sed "s:latsul:${south}:g"   \
  > ${Eta_ucl}/CTLTEMPLATE_FIXED
#
cat ${Eta_ucl}/cntrl.parm_NOPACK_temp | sed "s:ptx:$imout1:g"       \
                                      | sed "s:pty:$jmout1:g"       \
                                      | sed "s:reseta:${res}:g"     \
                                      | sed "s:lonoeste:${west2}:g" \
                                      | sed "s:latsul:${south}:g"   \
  > ${Eta_ucl}/cntrl.parm_NOPACK
