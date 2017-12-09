#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --time=0-00:20:00 # run time in days-hh:mm:ss
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --error=/srv/home/nroberts5/FinalProject/error_out2.err
#SBATCH --output=/srv/home/nroberts5/FinalProject/output2.out
#SBATCH --gres=gpu:0
module load cuda
./least_squares_test 1 10000
./least_squares_test 5 10000
./least_squares_test 10 10000
./least_squares_test 20 10000
./least_squares_test 30 10000
./least_squares_test 40 10000