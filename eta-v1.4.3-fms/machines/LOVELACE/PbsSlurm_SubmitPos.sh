#PBS -l nodes=NNODES:ppn=NCPUS
#PBS -l walltime=04:50:00
#PBS -q QFNAME
#PBS -N SPNAME
#PBS -o EONAME.out
#PBS -e EONAME.err
######################
module load openmpi/4.1.1-gcc-9.4.0
