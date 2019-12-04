//based on the example from https://www.codewithc.com/matrix-multiplication-in-c/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#define SEED 539147254

int** allocMatrix(int N, int M);
void fillMatrix(int **matrix, int N, int M);
void printMatrix(int **matrix, int N, int M);
void releaseMatrix(int **matrix);

int main(int argc, char *argv[]) {

	int column_first, row_first, column_second, row_second, sum = 0;
	int** first_array, **second_array, **result_array;

	if (argc < 5) {
		printf("Matrix1: first parameter #rows, second #columns\n"
				"Matrix2: third parameter #rows, fourth #columns\n");
		return EXIT_FAILURE;
	}

	row_first = atoi(argv[1]);
	column_first = atoi(argv[2]);
	row_second = atoi(argv[3]);
	column_second = atoi(argv[4]);

	//Checking if Matrix Multiplication is possible
	if (column_first != row_second) {
		printf(
				"Matrices with entered orders can't be multiplied with each other.\n"
						"The column of first matrix should be equal to row of second.\n");
		return EXIT_FAILURE;
	}

	first_array = allocMatrix(row_first, column_first);
	second_array = allocMatrix(row_second, column_second);
	result_array = allocMatrix(row_first, column_second);

	fillMatrix(first_array, row_first, column_first);
	fillMatrix(second_array, row_second, column_second);

	//printf("\nFirst matrix:\n");
	//printMatrix(first_array, row_first, column_first);

	//printf("\nSecond matrix:\n");
	//printMatrix(second_array, row_second, column_second);

	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);

	//Carrying out matrix multiplication operation
	#pragma omp parallel for schedule(dynamic) firstprivate(row_first, row_second, column_second, first_array, second_array, sum) shared(result_array) default(none)
	for (int i = 0; i < row_first; i++) {
		for (int j = 0; j < column_second; j++) {
			for (int k = 0; k < row_second; k++) {
				sum = sum + first_array[i][k] * second_array[k][j];
			}

			result_array[i][j] = sum;
			sum = 0;
		}
	}

	gettimeofday(&tv2, NULL);

    printf ("%f;\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));

	//Printing the final product matrix
	//printf("\nThe product of entered matrices is:\n");
	//printMatrix(result_array, row_first, column_second);

	releaseMatrix(result_array);
	releaseMatrix(first_array);
	releaseMatrix(second_array);

	return EXIT_SUCCESS;
}

int** allocMatrix(int N, int M) {
	   int *data = (int *)malloc(N*M*sizeof(int));
	    int** array= (int **)malloc(M*sizeof(int*));
	    for (int i=0; i<M; i++)
	        array[i] = &(data[N*i]);

	    return array;
}

void fillMatrix(int **matrix, int N, int M) {
	srand(SEED);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			matrix[i][j] = (int) rand() % 10;
		}
	}
}

void printMatrix(int **matrix, int N, int M) {

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			printf("%d\t", matrix[i][j]);
		}
		printf("\n");
	}
}

void releaseMatrix(int **matrix) {
	free(matrix[0]);
	free(matrix);
}
