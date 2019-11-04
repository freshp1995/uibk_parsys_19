CC=gcc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm
MPI=mpicc

.PHONEY: all
all: heat_stencil_3D_cubes
all: heat_stencil_3D_slabs


heat_stencil_3D_cubes: heat_stencil_3D_cubes.c
	@$(MPI) $(CC_FLAGS) $< -o heat_stencil_3D_cubes

heat_stencil_3D_slabs: heat_stencil_3D_slabs.c
	@$(MPI) $(CC_FLAGS) $< -o heat_stencil_3D_slabs


.PHONEY: clean
clean:
	@rm heat_stencil_3D_cubes
	@rm heat_stencil_3D_slabs


.PHONEY: run
run: all
	@echo "cubes 3D:"
	@./heat_stencil_3D_cubes
	@echo
	@echo "slabs 3D:"
	@./heat_stencil_3D_slabs
	@echo
