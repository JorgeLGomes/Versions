#!/bin/bash -x
get_gdrive.sh 1P8vxBwabdokJAcJoAFc9_1oOTg6CvQNM  ATMOS_CFS.01_Exp1.tgz
tar -zxvf ATMOS_CFS.01_Exp1.tgz
get_gdrive.sh 1x4FcNr2x6-gMldD6Du0VFF5j-UPMlYPC  SST_CFS.01_Exp1.tgz
tar -zxvf SST_CFS.01_Exp1.tgz
get_gdrive.sh 1MgI2rl6VCA0EfkgDwB2i-WG5LKavDjUW  ATMOS_CFS.01_Exp2.tgz
tar -zxvf ATMOS_CFS.01_Exp2.tgz
get_gdrive.sh 1x4FcNr2x6-gMldD6Du0VFF5j-UPMlYPC  SST_CFS.01_Exp2.tgz
tar -zxvf SST_CFS.01_Exp2.tgz
