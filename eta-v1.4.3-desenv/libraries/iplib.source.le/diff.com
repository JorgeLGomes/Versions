#! /bin/ksh -f

pwd

list=`ls *.f`
for name in $list
do
diff $name ../iplib.source/$name
echo "===== done $name ===== "
done
