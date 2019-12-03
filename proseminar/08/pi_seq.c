/* Program to compute Pi using Monte Carlo methods
 * based on example program from https://www.dartmouth.edu/~rc/classes/soft_dev/C_simple_ex.html
 * Authors: Raphael Gruber, Patrick Lanzinger
 * */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

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

	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);
	for (int i = 0; i < niter; i++) {
		x = (double) rand() / RAND_MAX;
		y = (double) rand() / RAND_MAX;
		z = x * x + y * y;
		if (z <= 1)
			count++;
	}
	pi = (double) count / niter * 4;
	//printf("# of samples= %d , estimate of pi is %g \n", niter, pi);
	//printf("%d; %g\n", niter, pi);

	gettimeofday(&tv2, NULL);

    printf ("%d;%f;\n",pi, (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));

	return EXIT_SUCCESS;
}
