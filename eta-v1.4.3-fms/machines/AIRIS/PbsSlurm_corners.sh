#SBATCH -o EXPPATH/Corners.o
#SBATCH -e EXPPATH/Corners.e
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH -J corners
#SBATCH -A inpe
#SBATCH -p processing
module load openmpi/4.1.2-intel
module load nvhpc/23.1
