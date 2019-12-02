/* Program to compute Pi using Monte Carlo methods
 * based on example program from https://www.dartmouth.edu/~rc/classes/soft_dev/C_simple_ex.html
 * Authors: Raphael Gruber, Patrick Lanzinger
 * */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define SEED 35791246

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("no number of samples given (program parameter)\n");
		return EXIT_FAILURE;
	}

	int niter = atoi(argv[1]);
	double x, y;
	long count = 0; /* # of points in the 1st quadrant of unit circle */
	double z;
	double pi;

	/* initialize random numbers */
	srand(SEED);
	count = 0;

	#pragma omp parallel for reduction(+:count) schedule(dynamic) shared(niter) private(x,y,z) default(none) 
	for (int i = 0; i < niter; i++) {
		x = (double) rand() / RAND_MAX;
		y = (double) rand() / RAND_MAX;
		z = x * x + y * y;
		if (z <= 1)
			count += 1;
	}
	pi = (double) count / niter * 4;
	//printf("# of samples= %d , estimate of pi is %g \n", niter, pi);
	printf("%d; %g\n", niter, pi);

	return EXIT_SUCCESS;
}
