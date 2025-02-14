#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

typedef double value_t;

// -- Matrix utilities --

typedef value_t ***Matrix;

typedef value_t **Vector;

Matrix createMatrix(int N);

Vector createVector(int N);

void releaseMatrix(Matrix m, int size);

void releaseVector(Vector m, int size);

void fill_Matrix(Matrix m, int size, int x, int y, int z);

int verify(Matrix m, int size);

void getCell(int* start,int* end, Matrix m, int size, Vector output);

void insertCell(value_t **c, int y, int z, Matrix m, int size);

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
		
	if(roundf(cbrt(numProcs)) != cbrt(numProcs)) {
		printf("Cubic root of rank size must be a natural number\n");
		MPI_Finalize();
		return EXIT_FAILURE;
	}


	int T = N_big * 500;
	
	//+4 because left and right behind before ghost cell
	int N = (N_big/cbrt(numProcs)) + 4;
	
	Matrix A = createMatrix(N);
	Matrix B = createMatrix(N);

	//fill the array
	int X = N_big / 4;
	int Y = N_big / 5;
	int Z = N_big / 2;
	
	MPI_Comm newComm;
	
    int dims[] = {cbrt(numProcs),cbrt(numProcs),cbrt(numProcs)};
    
    //int dims[3] = {0,0,0};

	//MPI_Dims_create(numProcs, 3, dims);
    
    int periods[] = {1,1,1};
    int reorder = 1;
    
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &newComm);
    
	//MPI_Comm_rank(newComm, &rank);
	
	
	if(rank == 0) {
		printf("Computing heat-distribution for romm size %dx%dx%d for %d timestamps\n", N_big, N_big, N_big, T);

		//create a buffer
		Matrix A_big = createMatrix(N_big);

	 
		fill_Matrix(A_big, N_big, X, Y, Z);

		//create a second buffer for the computation
		Matrix B_big = createMatrix(N_big);
		
		for (int i = 0; i < N_big; i++) {
			for (int j = 0; j < N_big; j++) {
            for (int k = 0; k < N_big; k++) {
				//printf("A_big%f --> x: %d y: %d z: %d\n", A_big[i][j][k], j, i, k);
		}}}
		
		for(int i = 0; i < numProcs; i++) {
			
			MPI_Datatype subArray;
			int coords[3];
			MPI_Cart_coords(newComm, i, 3, coords);
				
			int array_size[] = {N_big,N_big, N_big};
			int array_subsize[] = {N-4,N-4,N-4};
			int array_start[] = {coords[0]*array_subsize[0],coords[1]*array_subsize[1],coords[2]*array_subsize[2]}; 

			MPI_Type_create_subarray(3, array_size, array_subsize, array_start, MPI_ORDER_C, MPI_DOUBLE, &subArray);
			MPI_Type_commit(&subArray);
			
			MPI_Send(&(A_big[0][0][0]), 1, subArray, i, 0, newComm);
			MPI_Send(&(B_big[0][0][0]), 1, subArray, i, 0, newComm);
			
			MPI_Type_free(&subArray);
				
		}
		
		releaseMatrix(A_big, N_big);
		releaseMatrix(B_big,N_big);
    
	}

	//A[1][1][1] - first cell left out for ghost cells
	MPI_Recv(&(A[1][1][1]), (N-4)*(N-4)*(N-4), MPI_DOUBLE, 0, 0, newComm, MPI_STATUS_IGNORE);
    MPI_Recv(&(B[1][1][1]), (N-4)*(N-4)*(N-4), MPI_DOUBLE, 0, 0, newComm, MPI_STATUS_IGNORE);
    
    printf("Recive subarray size %d \n", (N-4)*(N-4)*(N-4));
    
	
/*
    value_t ghost_left[N][N];
    value_t ghost_right[N][N];
    value_t ghost_up[N][N];
    value_t ghost_down[N][N];
    value_t ghost_behind[N][N];
    value_t ghost_before[N][N];*/
    
    
    Vector ghost_left = createVector(N);
    Vector ghost_right = createVector(N);
    Vector ghost_up = createVector(N);
    Vector ghost_down = createVector(N);
    Vector ghost_behind = createVector(N);
    Vector ghost_before = createVector(N);
    
   
    
	int left_rank;
	int right_rank;
	int up_rank;
	int down_rank;
	int behind_rank;
	int before_rank;
	MPI_Cart_shift(newComm, 0, 1, &left_rank, &right_rank);
	MPI_Cart_shift(newComm, 1, 1, &up_rank, &down_rank);
	MPI_Cart_shift(newComm, 2, 1, &before_rank, &behind_rank);

	
    
    time_t start = clock();
    //for each time step
	for (int t = 0; t < T; t++) {
		MPI_Send(&A[1], N, MPI_DOUBLE, up_rank, 0, newComm);
		MPI_Send(&A[N-2], N, MPI_DOUBLE, down_rank, 0, newComm);
		
		int startLeft[3] = {1,1,1};
		int endLeft[3] = {2,N-2,N-2};		
		Vector tempArray = createVector(N);
		getCell(startLeft,endLeft,A,N-2, tempArray);
		MPI_Send(&tempArray, N-2, MPI_DOUBLE, left_rank, 0, newComm);

		int startRight[3] = {N-1,1,1};
		int endRight[3] = {N-2,N-2,N-2};	
		getCell(startRight,endRight,A,N-2,tempArray);
		MPI_Send(&tempArray, N-2, MPI_DOUBLE, right_rank, 0, newComm);
		
		int startBehind[3] = {1,1,N-2};
		int endBehind[3] = {N-2,N-2,N-1};
		getCell(startBehind,endBehind,A,N-2,tempArray);
		MPI_Send(&tempArray, N-2, MPI_DOUBLE, behind_rank, 0, newComm);
		
		int startBefore[3] = {1,1,1};
		int endBefore[3] = {N-2,N-2,2};
		getCell(startBefore,endBefore,A,N-2,tempArray);
		MPI_Send(&tempArray, N-2, MPI_DOUBLE, before_rank, 0, newComm);
		
		MPI_Recv(ghost_up, N, MPI_DOUBLE, up_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_down, N, MPI_DOUBLE, down_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_left, N, MPI_DOUBLE, left_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_right, N, MPI_DOUBLE, right_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_behind, N, MPI_DOUBLE, behind_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_before, N, MPI_DOUBLE, before_rank, 0, newComm, MPI_STATUS_IGNORE);

		printf("finished\n");

		A[0] = ghost_up;
		A[N-1] = ghost_down;
		B[0] = ghost_up;
		B[N-1] = ghost_down;
		insertCell(ghost_left, 0,0, A, N);
		insertCell(ghost_left, 0,0, B, N);
		insertCell(ghost_right, N,0, A, N);
		insertCell(ghost_right, N,0, B, N);
		insertCell(ghost_behind, 0,N, A, N);
		insertCell(ghost_behind, 0,N, B, N);
		insertCell(ghost_before, N,N, A, N);
		insertCell(ghost_before, N,N, B, N);
		printf("done\n");


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
                    value_t t_above = (i != 0) ? A[i - 1][j][k] : tc;
					value_t t_left = (j != 0) ? A[i][j - 1][k] : tc;
					value_t t_right = (j != N - 1) ? A[i][j + 1][k] : tc;
					value_t t_below = (i != N - 1) ? A[i + 1][j][k] : tc;


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
			if (!(t % 1000)) {
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
            y[i][j] = malloc(sizeof(value_t) * N);

            if (y[i][j] == NULL) {
                perror("Could not allocate memory");
                return NULL;
            }
        }
    }

    return y;
}

void releaseMatrix(Matrix m, int size) {

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            free(m[i][j]);
        }
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

    m[y][x][z] = 273.0 + 60.0;

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



void insertCell(value_t **c, int y, int z, Matrix m, int size) {
	
	
	for(int i = size-2; i >= 0; i--) {
		for(int j = size-1; j >= 0; j--) {
			for(int k = size-1; k >=0; k--) {
				m[i][j+1][k+1] = m[i][j][k];
				if(j == y && k == z) {
					m[i][j][k] = c[j][k];
					break;
				}
			}
		}
	
	/*for(int i = size-2; i >= 0; i--) {
		for(int j = size-1; j >= 0; j--) {
			for(int k = size-1; k >=0; k--) {
				m[i][j+1][k+1] = m[i][j][k];
				if(j == y && k == z) {
					m[i][j][k] = c[j][k];
					break;
				}
			}
		}
	}*/
	
	
}

Vector createVector(int N) {
    Vector y = malloc(sizeof(value_t) * N);

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
    }

    return y;
}

void releaseVector(Vector m, int size) {
    for (int i = 0; i < size; i++) {
        free(m[i]);
    }

    free(m);
}
