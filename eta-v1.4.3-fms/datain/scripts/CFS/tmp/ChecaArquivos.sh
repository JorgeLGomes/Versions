#!/bin/bash -x

foreach i in $(ls -l /dados/grpeta/dsk003/dados/CFS/output/2022011512_M01/CFS_2022011512E.????); do
echo $i |awk "(print $5)"
done
