#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --time=0-00:20:00 # run time in days-hh:mm:ss
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --error=/srv/home/nroberts5/FinalProject/error_out.err
#SBATCH --output=/srv/home/nroberts5/FinalProject/output.out
#SBATCH --gres=gpu:0
module load cuda
./NLSQ 40 10000 1
