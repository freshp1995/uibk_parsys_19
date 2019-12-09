CC=gcc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm -fopenmp

.PHONEY: all
all: n-queen
all: n-queen_opt
all: n-queen_mp


n-queen: n-queen.c
	@$(CC) $(CC_FLAGS) $< -o n-queen

n-queen_opt: n-queen_opt.c
	@$(CC) $(CC_FLAGS) $< -o n-queen_opt
	
n-queen_mp: n-queen_mp.c
	@$(CC) $(CC_FLAGS) $< -o n-queen_mp





.PHONEY: clean
clean:
	@rm n-queen
	@rm n-queen_opt
	@rm n-queen_mp



.PHONEY: run
run: all
	@echo "queen"
	@./n-queen
	@echo
	@echo "queen_opt"
	@./n-queen_opt
	@echo
	@echo "queen_mp"
	@./n-queen_mp
	@echo
