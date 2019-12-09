CC=gcc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm -fopenmp

.PHONEY: all
all: n-queen
all: n-queen_opt
all: n-queen_mp


queen: queen.c
	@$(CC) $(CC_FLAGS) $< -o queen

queen_opt: queen_opt.c
	@$(CC) $(CC_FLAGS) $< -o queen_opt
	
queen_mp: queen_mp.c
	@$(CC) $(CC_FLAGS) $< -o queen_mp





.PHONEY: clean
clean:
	@rm queen
	@rm queen_opt
	@rm queen_mp



.PHONEY: run
run: all
	@echo "queen"
	@./queen
	@echo
	@echo "queen_opt"
	@./queen_opt
	@echo
	@echo "queen_mp"
	@./queen_mp
	@echo
