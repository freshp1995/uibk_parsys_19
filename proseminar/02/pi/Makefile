CC=gcc
MPI=mpicc
CC_FLAGS=-O2 -std=c99 -Wall -Wextra -pedantic

.PHONEY: all
all: pi_seq
all: pi_mpi

pi_seq: pi_seq.c
	@$(CC) $(CC_FLAGS) $< -o pi_seq

pi_mpi: pi_seq.c
	@$(MPI) $(CC_FLAGS) $< -o pi_mpi

.PHONEY: clean
clean:
	@rm pi_seq
	@rm pi_mpi

