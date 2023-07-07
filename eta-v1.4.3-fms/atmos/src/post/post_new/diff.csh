#!/bin/csh -x
mkdir ./DIFF
foreach i (*.f) 
sdiff -b -B -s $i /lustre_xc50/d_rodrigues/versoes/eta-v1.4.3-desenv/atmos/src/post/post_new/$i > ./DIFF/$i.diff
end
