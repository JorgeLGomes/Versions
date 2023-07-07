#SBATCH -o ${Eta_run}/Preproc.o
#SBATCH -e ${Eta_run}/Preproc.e
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH -J ${Exp}P${Run_Date}
#SBATCH -A inpe
#SBATCH -p processing
module load openmpi/4.1.2-intel
module load nvhpc/23.1

