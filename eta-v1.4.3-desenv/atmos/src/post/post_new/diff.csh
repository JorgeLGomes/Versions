#!/bin/csh -x
mkdir ./DIFF
foreach i (*.f) 
sdiff -b -B -s $i /lustre_xc50/gustavo_sueiro/Oper/worketa/eta/src/post/post_new/$i > ./DIFF/$i.diff
end
