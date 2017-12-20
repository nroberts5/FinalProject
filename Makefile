# Warnings
WFLAGS	:= -Wall -Wextra -Wsign-conversion -Wsign-compare

# Optimization and architecture
OPT		:= -O3
ARCH   	:= -march=native

# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++11

# Linker options
LDOPT 	:= $(OPT)
LDFLAGS := 
BIN = "/usr/local/gcc/6.4.0/bin/gcc"
.DEFAULT_GOAL := all

.PHONY: debug
debug : OPT  := -O0 -g -fno-omit-frame-pointer -fsanitize=address
debug : LDFLAGS := -fsanitize=address
debug : ARCH :=
debug : $(EXEC)

all : MC_NLSQ MC_NLSQ_MPI

MC_NLSQ : MC_NLSQ.cu
	@ module load cuda;nvcc -o MC_NLSQ $(OPT) -Xcompiler -fopenmp MC_NLSQ.cu -ccbin $(BIN)

MC_NLSQ_MPI: MC_NLSQ_MPI.cpp
	@ module load openmpi/2.1.1;mpicxx -o MC_NLSQ_MPI $(CXXSTD) $(OPT) -fopenmp  MC_NLSQ_MPI.cpp


.PHONY: clean
clean:
	rm -f MC_NLSQ MC_NLSQ_MPI
	rm -f *.err *.out
