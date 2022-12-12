#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=60G
#SBATCH --partition=gpu
#SBATCH --qos=gpu_free
#SBATCH --account=cadmos
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

module purge
module load gcc cuda matlab

srun -n 1 matlab -nodisplay -nosplash -nodesktop -r main_HCP_CUDA
