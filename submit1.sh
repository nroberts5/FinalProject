#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --time=0-00:05:00 # run time in days-hh:mm:ss
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --error=/srv/home/nroberts5/759--nroberts5/HW10/error_out.err
#SBATCH --output=/srv/home/nroberts5/759--nroberts5/HW10/output.out
#SBATCH --gres=gpu:1
module load cuda
./problem1 $1
# ./problem1 2
# ./problem1 4
# ./problem1 8
# ./problem1 16
# ./problem1 32
# ./problem1 64
# ./problem1 128
# ./problem1 256
# ./problem1 512
# ./problem1 1024
# ./problem1 2048
# ./problem1 4096
# ./problem1 8192
# ./problem1 16384
# ./problem1 32768
# ./problem1 65536
# ./problem1 131072
# ./problem1 262144
# ./problem1 524288
# ./problem1 1048576
# ./problem1 2097152
# ./problem1 4194304
# ./problem1 8388608
# ./problem1 16777216