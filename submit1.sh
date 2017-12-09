#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --time=0-00:20:00 # run time in days-hh:mm:ss
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --error=/srv/home/nroberts5/FinalProject/error_out1.err
#SBATCH --output=/srv/home/nroberts5/FinalProject/output1.out
#SBATCH --gres=gpu:0
module load cuda
./least_squares_test 1 10
./least_squares_test 5 10
./least_squares_test 10 10
./least_squares_test 20 10
./least_squares_test 30 10
./least_squares_test 40 10

./least_squares_test 1 100
./least_squares_test 5 100
./least_squares_test 10 100
./least_squares_test 20 100
./least_squares_test 30 100
./least_squares_test 40 100

./least_squares_test 1 1000
./least_squares_test 5 1000
./least_squares_test 10 1000
./least_squares_test 20 1000
./least_squares_test 30 1000
./least_squares_test 40 1000
