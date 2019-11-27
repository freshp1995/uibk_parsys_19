CC=gcc
CC_FLAGS=-O0 -g -ggdb -std=c99 -Wall -Wextra -pedantic -lm -fopenmp

.PHONEY: all
all: n-body-seq
all: n-body
all: n-body-2


n-body-seq: n-body-seq.c
	@$(CC) $(CC_FLAGS) $< -o n-body-seq

n-body: n-body.c
	@$(CC) $(CC_FLAGS) $< -o n-body

n-body-2: n-body-2.c
	@$(CC) $(CC_FLAGS) $< -o n-body-2



.PHONEY: clean
clean:
	@rm n-body-seq
	@rm n-body
	@rm n-body-2


.PHONEY: run
run: all
	@echo "n-body-seq"
	@./n-body-seq
	@echo
	@echo "n-body"
	@./n-body
	@echo
	@echo "n-body-2"
	@./n-body-2
	@echo