//from https://rosettacode.org/wiki/N-queens_problem#C

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
 
int count = 0;
int solve(int n, int col, int *hist, int print)
{
	if (col == n) {
		count++;
		if(print) {
			printf("\nNo. %d\n-----\n", count);
			for (int i = 0; i < n; i++, putchar('\n'))
				for (int j = 0; j < n; j++)
					putchar(j == hist[i] ? 'Q' : ((i + j) & 1) ? ' ' : '.');
 
		}
		return count;
	}
 
#define attack(i, j) (hist[j] == i || abs(hist[j] - i) == col - j)
	for (int i = 0, j = 0; i < n; i++) {
		for (j = 0; j < col && !attack(i, j); j++);
		if (j < col) 
			continue;
 
		hist[col] = i;
		solve(n, col + 1, hist, print);
	}
	
	return count;
	
}
 
int main(int argc, char **argv)
{

	
	if(argc < 2) {
    printf("usage: number of queens, optionally 1 for print out\n");
    return 1;
	}
  
	int n = atoi(argv[1]);

	int print = 0;
	if (argc == 3) 
		print = atoi(argv[2]);


	
	int hist[n];
	
	
	clock_t begin = clock();
  	int result = solve(n, 0, hist,print);
  	clock_t end = clock();

  	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  	printf("%f\n", time_spent);

  	//printf("Number of solution for %d is %d\n",n,result);
	
	return EXIT_SUCCESS;

}
