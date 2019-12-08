CC=gcc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm -fopenmp

.PHONEY: all
all: n-queen
all: n-queen_opt


queen: queen.c
	@$(CC) $(CC_FLAGS) $< -o queen

queen_opt: queen_opt.c
	@$(CC) $(CC_FLAGS) $< -o queen_opt





.PHONEY: clean
clean:
	@rm queen
	@rm queen_opt



.PHONEY: run
run: all
	@echo "queen"
	@./queen
	@echo
	@echo "queen_opt"
	@./queen_opt
	@echo

	@echo
