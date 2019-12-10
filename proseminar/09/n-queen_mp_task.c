//from https://rosettacode.org/wiki/N-queens_problem#C

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
 
int count = 0;
int solve(int n, int col, int *hist)
{
	if (col == n) {
		
		#pragma omp atomic
		count++;

		return count;
	}
 
	#define attack(i, j) (hist[j] == i || abs(hist[j] - i) == col - j)

	for (int i = 0, j = 0; i < n; i++) {
		for (j = 0; j < col && !attack(i, j); j++);
		if (j < col) 
			continue;
 
		hist[col] = i;

		#pragma omp task shared(hist, count) firstprivate(col, n) default(none)
		solve(n, col + 1, hist);
		#pragma omp taskwait
	}

	return count;
	
}
 
int main(int argc, char **argv)
{

	
	if(argc < 2) {
    printf("usage: number of queens\n");
    return 1;
	}
  
	int n = atoi(argv[1]);
	
	int hist[n];
	
	int result = 0;
	#pragma omp parallel 
	{
		#pragma omp task
		result = solve(n, 0, hist);
		
	}
  	printf("Number of solution for %d is %d\n",n, result);
	
	return EXIT_SUCCESS;

}
