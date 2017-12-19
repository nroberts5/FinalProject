#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --time=0-00:20:00 # run time in days-hh:mm:ss
#SBATCH --nodes=2
#SBATCH --cpus-per-task=10
#SBATCH --error=/srv/home/nroberts5/FinalProject/error_out.err
#SBATCH --output=/srv/home/nroberts5/FinalProject/output.out
#SBATCH --gres=gpu:0

module load openmpi/2.1.1
mpiexec -np 2 ./MPI_Test
