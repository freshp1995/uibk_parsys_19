CC=gcc
CC_FLAGS=-O2 -std=c99 -Wall -Wextra -pedantic
MPI=mpicc

.PHONEY: all
all: heat_stencil_1D_seq
all: heat_stencil_1D_mpi

heat_stencil_1D_seq: heat_stencil_1D_seq.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_1D_seq

heat_stencil_1D_mpi: heat_stencil_1D_mpi.c
	@$(MPI) $(CC_FLAGS) $< -o heat_stencil_1D_mpi

.PHONEY: clean
clean:
	@rm heat_stencil_1D_seq
	@rm heat_stencil_1D_mpi

.PHONEY: run
run: all
	@echo "Sequential:"
	@./heat_stencil_1D_seq
	@echo
	@echo "MPI:"
	@./heat_stencil_1D_mpi
	@echo
