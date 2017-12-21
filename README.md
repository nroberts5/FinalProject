# NRoberts_ECE759_FinalProject

## Without MPI
* Usage: ./MC_NLSQ NUMTHREADS Number_of_Sims

##With MPI
* Usage: mpiexec -np 2 ./MC_NLSQ_MPI NUMTHREADS Number_of_Sims

## Definitions:
* NUMTHREADS = Integer Number of OpenMP Threads Per Node
* Number_of_Sims = Integer Number of Monte Carlo Simulations to Perform

## Returns:
* 8 Files with simulated estimates (.out files)
* To standard output: execution time details.

## Examples of how I tested these codes (which worked for me) are in the files:

* submit.sh (for MC_NLSQ)
* submitMPI.sh (for MC_NLSQ_MPI)
