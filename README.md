# FinalProject
NRoberts_ECE759_FinalProject
Usage: ./MC_NLSQ NUMTHREADS Number_of_Sims
Usage: mpiexec -np 2 ./MC_NLSQ_MPI NUMTHREADS Number_of_Sims

Definitions:
NUMTHREADS = Integer Number of OpenMP Threads Per Node
Number_of_Sims = Integer Number of Monte Carlo Simulations to Perform

Returns:
8 Files with simulated estimates (.out files)

To standard output: execution time details.
