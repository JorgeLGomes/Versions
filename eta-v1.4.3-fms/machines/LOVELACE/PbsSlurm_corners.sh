#PBS -l nodes=1:ncpus=1
#PBS -l walltime=0:10:00
#PBS -q serial
#PBS -N  corners
#PBS -o EXPPATH/Corners.o
#PBS -e EXPPATH/Corners.e
module load openmpi/4.1.1-gcc-9.4.0
