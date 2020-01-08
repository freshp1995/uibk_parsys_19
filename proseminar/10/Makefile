CC=gcc
MPI=mpicc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm -fopenmp

.PHONEY: all
all: heat_stencil_2D


heat_stencil_2D: heat_stencil_2D.c
	@$(MPI) $(CC_FLAGS) $< -o heat_stencil_2D



.PHONEY: clean
clean:
	@rm heat_stencil_2D



.PHONEY: run
run: all
	@echo "queen"
	@./heat_stencil_2D
	@echo
