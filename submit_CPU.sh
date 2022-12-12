#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mem=60G
#SBATCH --account=cadmos
#SBATCH -p debug
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

module load matlab

srun -n 1 matlab -nodisplay -nosplash -nodesktop -r main_HCP
