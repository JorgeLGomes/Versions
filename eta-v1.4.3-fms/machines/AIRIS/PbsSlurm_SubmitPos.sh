#SBATCH -o ${Eta_run}/EONAME.o
#SBATCH -e ${Eta_run}/EONAME.e
#SBATCH --ntasks-per-node=NCPUS
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH -J SPNAME
#SBATCH -A inpe
#SBATCH -p processing
module load openmpi/4.1.2-intel
module load nvhpc/23.1
