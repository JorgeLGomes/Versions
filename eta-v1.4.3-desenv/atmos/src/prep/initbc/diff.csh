#!/bin/csh -x

foreach i (*.f) 
sdiff -b -B -s $i /lustre_xc50/gustavo_sueiro/Oper/worketa/eta/src/prep/initbc/$i > ./DIFF/$i.diff
end
