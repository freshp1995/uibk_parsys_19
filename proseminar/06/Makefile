CC=gcc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm
MPI=mpicc

.PHONEY: all
all: n-body



n-body: n-body.c
	@$(CC) $(CC_FLAGS) $< -o n-body



.PHONEY: clean
clean:
	@rm n-body


.PHONEY: run
run: all
	@echo "n-body"
	@./n-body
	@echo