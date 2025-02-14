/* Program to compute Pi using Monte Carlo methods
 * based on example program from https://www.dartmouth.edu/~rc/classes/soft_dev/C_simple_ex.html
 * Authors: Raphael Gruber, Patrick Lanzinger
 * */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#define SEED 35791246

int main(int argc, char **argv) {
	if (argc < 2) {
		printf("no number of samples (program parameter)");
		return EXIT_FAILURE;
	}

	int numProcs, rank;

	/* initialize mpi*/
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

	int niter = atoi(argv[1]);
	double x, y;
	long count = 0; /* # of points in the 1st quadrant of unit circle */
	double z;
	double pi = 0;

	/* initialize random numbers */
	srand(SEED + rank);
	count = 0;



	int from = rank * niter / numProcs;
	int to = (rank + 1) * niter / numProcs;

	for (int i = from; i < to; i++) {
		x = (double) rand() / RAND_MAX;
		y = (double) rand() / RAND_MAX;
		z = x * x + y * y;
		if (z <= 1)
			count++;
	}

	long allCount = 0;

	/* sum up all counts */
	MPI_Reduce(&count, &allCount, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		pi = (double) allCount / niter * 4;

		//printf("# of samples= %d , estimate of pi is %g \n", niter, pi);
		printf("%d; %g\n", niter, pi);

	}

	MPI_Finalize();

	return EXIT_SUCCESS;

}
