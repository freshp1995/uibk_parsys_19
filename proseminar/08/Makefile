CC=gcc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm -fopenmp

.PHONEY: all
all: pi_seq
all: pi_mp
all: matrix_seq
all: matrix_mp


pi_seq: pi_seq.c
	@$(CC) $(CC_FLAGS) $< -o pi_seq

pi_mp: pi_mp.c
	@$(CC) $(CC_FLAGS) $< -o pi_mp

matrix_seq: matrix_seq.c
	@$(CC) $(CC_FLAGS) $< -o matrix_seq

matrix_mp: matrix_mp.c
	@$(CC) $(CC_FLAGS) $< -o matrix_mp




.PHONEY: clean
clean:
	@rm pi_seq
	@rm pi_mp
	@rm matrix_seq
	@rm matrix_mp


.PHONEY: run
run: all
	@echo "pi_seq"
	@./pi_seq
	@echo
	@echo "pi_mp"
	@./pi_mp
	@echo
	@echo "matrix_seq"
	@./matrix_seq
	@echo
	@echo "matrix_mp"
	@./matrix_mp
	@echo
