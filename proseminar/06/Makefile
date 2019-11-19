CC=gcc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm
MPI=mpicc

.PHONEY: all
all: n-body



n-body: n-body.c
	@$(MPI) $(CC_FLAGS) $< -o exec



.PHONEY: clean
clean:
	@rm exec


.PHONEY: run
run: all
	@echo "n-body"
	@./n-body
	@echo