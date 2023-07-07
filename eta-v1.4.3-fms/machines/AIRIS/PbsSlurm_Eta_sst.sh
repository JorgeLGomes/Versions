#SBATCH -o ${Eta_run}/SST.o
#SBATCH -e ${Eta_run}/SST.e
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH -J EtaSST
#SBATCH -A inpe
#SBATCH -p processing
module load openmpi/4.1.2-intel
module load nvhpc/23.1
