#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --time=0-00:05:00 # run time in days-hh:mm:ss
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --error=/srv/home/nroberts5/FinalProject/error_out.err
#SBATCH --output=/srv/home/nroberts5/FinalProject/output.out
#SBATCH --gres=gpu:0
module load cuda
./main 20 1000
./main 20 10000
./main 20 100000
# ./main 20 1000000
# ./main 20 10000000