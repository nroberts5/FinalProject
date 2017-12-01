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

all : main

# problem1: problem1.cu 
# 	module load cuda;nvcc -o problem1 $(OPT) problem1.cu -ccbin $(BIN)

main : main.cu
	@ module load cuda;nvcc -o main $(OPT) main.cu -ccbin $(BIN)

# TODO: add targets for building executables

.PHONY: clean
clean:
	rm -f main
	rm -f *.err *.out
