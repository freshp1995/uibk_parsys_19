#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

typedef double value_t;

// -- Matrix utilities --

typedef value_t ***Matrix;

typedef value_t **Vector;

Matrix createSlabMatrix(int height, int N);

Matrix createMatrix(int N);

Vector createVector(int N);

void releaseMatrix(Matrix m, int size);

void releaseVector(Vector m, int size);

void fill_Matrix(Matrix m, int size, int x, int y, int z);

int verify(Matrix m, int size);

void getCell(int* start,int* end, Matrix m, int size, Vector output);

void insertCell(value_t **c, int y, int z, Matrix m, int size);

void printVector (Vector m, int size, int rank);
	
void printMatrix (Matrix m, int size);

// -- simulation code ---
long timediff(clock_t t1, clock_t t2) {
    long elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000;
    return elapsed;
}

int main(int argc, char **argv) {
    int N_big = 10;
    if (argc > 1) {
        N_big = atoi(argv[1]);
    }
    
    int rank, numProcs;
        
    /* initialize mpi*/
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		
	if(roundf(N_big/numProcs) != (N_big/numProcs)) {
		printf("Cubic root of rank size must be a natural number\n");
		MPI_Finalize();
		return EXIT_FAILURE;
	}


	int T = N_big * 500;
	
	int N = (N_big/numProcs);
	
	Matrix A = createSlabMatrix(N, N_big);
	Matrix B = createSlabMatrix(N, N_big);


	//fill the array
	int X = N_big / 4;
	int Y = N_big / 5;
	int Z = N_big / 2;

	
	
	if(rank == 0) {
		printf("Computing heat-distribution for romm size %dx%dx%d for %d timestamps\n", N_big, N_big, N_big, T);

		//create a buffer
		Matrix A_big = createMatrix(N_big);

	 
		fill_Matrix(A_big, N_big, X, Y, Z);
		
		
		for(int i = 1; i < numProcs; i++) {
			
	
			MPI_Send(&(A_big[i*N][0][0]), N_big*N_big*N, MPI_DOUBLE, i, 42, MPI_COMM_WORLD);
				
		}

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N_big; j++) {
                for (int k = 0; k < N_big; k++) {
                    A[i][j][k] = A_big[i][j][k];
                }
            }
        }

		releaseMatrix(A_big, N_big);    
	}

    

    if (rank != 0) {
        printf("Ready to receive data %d\n", rank);
	    MPI_Recv(&(A[0][0][0]), N_big*N_big*N, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Recive subarray size %d \n", (N_big)*(N_big)*(N));
    }

    
    	
	MPI_Comm newComm;
	
    int dims[] = {numProcs};
    
    int periods[] = {1};
    int reorder = 0;
    
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &newComm);
    
	//MPI_Comm_rank(newComm, &rank);
 
    
    Vector ghost_up = createVector(N_big);
    Vector ghost_down = createVector(N_big);
    
   
    

	int up_rank;
	int down_rank;
	MPI_Cart_shift(newComm, 0, 1, &up_rank, &down_rank);


    time_t start = clock();
    //for each time step
    printf("Begin to calculate--------");
	for (int t = 0; t < T; t++) {
		

        if (rank % 2 == 0) {
            MPI_Send(&(A[0][0][0]), N_big*N_big, MPI_DOUBLE, up_rank, 1, newComm);
            MPI_Send(&(A[N-1][0][0]), N_big*N_big, MPI_DOUBLE, down_rank, 0, newComm);
            
            
            MPI_Recv(&(ghost_up[0][0]), N_big*N_big, MPI_DOUBLE, up_rank, 0, newComm, MPI_STATUS_IGNORE);
            MPI_Recv(&(ghost_down[0][0]), N_big*N_big, MPI_DOUBLE, down_rank, 1, newComm, MPI_STATUS_IGNORE);
        } else {

            MPI_Recv(&(ghost_down[0][0]), N_big*N_big, MPI_DOUBLE, down_rank, 1, newComm, MPI_STATUS_IGNORE);
            MPI_Recv(&(ghost_up[0][0]), N_big*N_big, MPI_DOUBLE, up_rank, 0, newComm, MPI_STATUS_IGNORE);
            

            MPI_Send(&(A[N-1][0][0]), N_big*N_big, MPI_DOUBLE, down_rank, 0, newComm);
            MPI_Send(&(A[0][0][0]), N_big*N_big, MPI_DOUBLE, up_rank, 1, newComm);
           
        }
	




        //we propagate the temparature
        for (long long i = 0; i < N; i++) {
            for (long long j = 0; j < N; j++) {
                for (long long k = 0; k < N; k++) {
                    if (i == Y && j == X && k == Z) {
                        B[i][j][k] = A[i][j][k];
                        continue;
                    }

                    //get temperatur at current position
                    value_t tc = A[i][j][k];


                    //get temperatur of adjacent cells
                    value_t t_behind = (k != 0) ? A[i][j][k - 1] : tc;
                    value_t t_before = (k != N - 1) ? A[i][j][k + 1] : tc;
                    value_t t_above = (i != 0) ? A[i - 1][j][k] : ghost_up[j][k];
					value_t t_left = (j != 0) ? A[i][j - 1][k] : tc;
					value_t t_right = (j != N - 1) ? A[i][j + 1][k] : tc;
					value_t t_below = (i != N - 1) ? A[i + 1][j][k] : ghost_down[j][k];


                    B[i][j][k] = tc + 0.1 * (
                                t_above + t_left + 
                                t_right + t_below + 
                                t_behind + t_before + ( -6 * tc )
                        );
                    //printf("%f --> x: %d y: %d z: %d\n", B[i][j][k] , j, i, k);
                }
            }
        }

        Matrix H = A;
        A = B;
        B = H;

		if(rank == 0) {
			if (!(t % 10)) {
				printf("Current timestamp t=%d\n", t);
			}
		}
    }
    time_t stop = clock();

    releaseMatrix(B, N);

    printf("Verification: %s\n", (verify(A, N)) ? "OK" : "FAILED");

	if(rank == 0) {
		long elapsed = timediff(start, stop);
		printf("elapsed: %ld ms\n", elapsed);
	}

    //release the Matrix again
    releaseMatrix(A, N);
    
     MPI_Finalize();
}

Matrix createMatrix(int N) {
    value_t *data = (value_t *)malloc(sizeof(value_t) * N * N * N);
    if (data == NULL) {
        perror("Could not allocate memory");
        return NULL;
    }

    Matrix y = malloc(sizeof(value_t) * N);

    if (y == NULL) {
        perror("Could not allocate memory");
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        y[i] = malloc(sizeof(value_t) * N);

        if (y[i] == NULL) {
            perror("Could not allocate memory");
            return NULL;
        }

        for (int j = 0; j < N; j++) {
            y[i][j] = &(data[N * (i * (int)pow(N, 1) + j * (int)pow(N, 0))]);
        }
    }

    return y;
}

Matrix createSlabMatrix(int height, int N) {
    value_t *data = (value_t *)malloc(sizeof(value_t) * height * N * N);
    if (data == NULL) {
        perror("Could not allocate memory");
        return NULL;
    }

    Matrix y = malloc(sizeof(value_t) * height);

    if (y == NULL) {
        perror("Could not allocate memory");
        return NULL;
    }

    for (int i = 0; i < height; i++) {
        y[i] = malloc(sizeof(value_t) * N);

        if (y[i] == NULL) {
            perror("Could not allocate memory");
            return NULL;
        }

        for (int j = 0; j < N; j++) {
            y[i][j] = &(data[N * (i * (int)pow(N, 1) + j * (int)pow(height, 0))]);
        }
    }

    return y;
}

void releaseMatrix(Matrix m, int size) {

    free(m[0][0]);

    for (int i = 0; i < size; i++) {
        free(m[i]);
    }

    free(m);
}

void fill_Matrix(Matrix m, int size, int x, int y, int z) {
    if (x > size || y > size) {
        return;
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                m[i][j][k] = 273.0;
            }
        }
    }

    m[0][0][0] = 273.0 + 60.0;

}

int verify(Matrix m, int size) {
    int success = 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                value_t temp = m[i][j][k];
                if (273 <= temp && temp <= 273 + 60) {
                    continue;
                }
                //printf("%f --> x: %d y: %d z: %d\n", temp, j, i, k);
                success = 0;
                break;
            }
        }
    }
    return success;
}

void getCell(int* start, int* end, Matrix m, int size, Vector output) {

	int tempX = 0, tempY = 0;
	
	for(int i = start[0]; i < end[0]; i++) {
		for(int j = start[1]; j < end[1]; j++) {
			for(int k = start[2]; k < end[2]; k++) {
				output[tempX][tempY] = m[i][j][k];
				tempX++;
				if(tempX >= size) {
					tempX = 0;
					tempY++;
				}
			}
		}
	}
	
	
}




Vector createVector(int N) {
   value_t *data = (value_t *)malloc(N*N*sizeof(value_t));
    Vector array= (int **)malloc(N*sizeof(value_t*));
    for (int i=0; i<N; i++)
        array[i] = &(data[N*i]);

    return array;
}

void releaseVector(Vector m, int size) {
free(m[0]);
free(m);
}


void printVector (Vector m, int size, int rank) {
	
	 for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
			printf("%f --> x: %d y: %d rank: %d\n", m[i][j], j, i, rank);
		}
	} printf("%d rank finished\n", rank);
	
}

void printMatrix (Matrix m, int size) {
	
	 for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
				printf("%f --> x: %d y: %d z: %d\n", m[i][j][k], j, i, k);
			}
		}
	}
	
}

void storeInBuffer(Vector v, value_t* b, int size) {
	int n = 0;
	for (int i = 0; i < size; i++) {
		for(int j = 0; j < size; j++) {
			b[n] = v[i][j];
			n++;
		}
  }
}
