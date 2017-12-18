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

all : least_squares_test MPI_Test

# problem1: problem1.cu 
# 	module load cuda;nvcc -o problem1 $(OPT) problem1.cu -ccbin $(BIN)

# main : main.cu
# 	@ module load cuda;nvcc -o main $(OPT) main.cu -ccbin $(BIN)

least_squares_test : least_squares_test.cu
	@ module load cuda;nvcc -o least_squares_test $(OPT) -Xcompiler -fopenmp least_squares_test.cu -ccbin $(BIN)

MPI_Test: MPI_Test.cpp
	@ module load openmpi/2.1.1;mpicxx -o MPI_Test $(CXXSTD) $(WFLAGS) $(OPT) -fopenmp  MPI_Test.cpp

# TODO: add targets for building executables

.PHONY: clean
clean:
	rm -f main least_squares_test MPI_Test
	rm -f *.err *.out
