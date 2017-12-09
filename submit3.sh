#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --time=0-00:20:00 # run time in days-hh:mm:ss
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --error=/srv/home/nroberts5/FinalProject/error_out3.err
#SBATCH --output=/srv/home/nroberts5/FinalProject/output3.out
#SBATCH --gres=gpu:0
module load cuda
./least_squares_test 40 100000 1
#./least_squares_test 1 100000
