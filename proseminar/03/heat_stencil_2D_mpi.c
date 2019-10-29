#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

typedef double value_t;

// -- vector utilities --

typedef value_t **Vector;

Vector createVector(int N);

void releaseVector(Vector m, int size);

void fill_vector(Vector m, int size, int x, int y);

int verify(Vector m, int size);

void printVector(Vector m, int size);

value_t *getColumn(int pos, Vector m, int size);

void insertColumn(value_t *c, int pos, Vector m, int size);

// -- simulation code ---

int main(int argc, char **argv) {
    int N_big = 10; 
    if (argc > 1) {
        N_big = atoi(argv[1]);
    }
    
    
	int rank;
	int numProcs;

	/* initialize mpi*/
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	
	if(numProcs%2 != 0) {
		printf("Rank size only multiple of 2");
		return EXIT_FAILURE;
	}
  
  

    int T = N_big * 500;
   	
	//+2 because left and right ghost cell
	int N = (N_big/sqrt(numProcs)) + 2;
	
	Vector A = createVector(N);
	Vector B = createVector(N);

	//fill the array
	int X = N_big / 4;
	int Y = N_big / 5;
	
	MPI_Comm newComm;
    int dims[] = {sqrt(numProcs),sqrt(numProcs)};
    int periods[] = {1,1};
    int reorder = 0;
	
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &newComm);

	if(rank == 0) {
		 printf("Computing heat-distribution for romm size %dx%d for %d timestamps\n", N_big, N_big, T);


		//create a buffer
		Vector A_big = createVector(N_big);


		fill_vector(A_big, N_big, X, Y);

		

		//create a second buffer for the computation
		Vector B_big = createVector(N_big);
		

		
		for(int i = 0; i < numProcs; i++) {
			MPI_Datatype subArray;
			int coords[2];
			MPI_Cart_coords(newComm, i, 2, coords);
			
			int array_size[] = {N_big,N_big};
			int array_subsize[] = {N-2,N-2};
			int array_start[] = {coords[0]*array_subsize[0],coords[1]*array_subsize[1]}; 

			MPI_Type_create_subarray(2, array_size, array_subsize, array_start, MPI_ORDER_C, MPI_DOUBLE, &subArray);
			MPI_Type_commit(&subArray);
			
			MPI_Send(&(A_big[0][0]), 1, subArray, i, 0, newComm);
			MPI_Send(&(B_big[0][0]), 1, subArray, i, 0, newComm);
			
		    MPI_Type_free(&subArray);
		    
		}
		
		
		
		releaseVector(A_big, N_big);
		releaseVector(B_big,N_big);
		

    
	}
	
	
	
	//A[1][1] - first cell left out for ghost cells
	MPI_Recv(&(A[1][1]), (N-2)*(N-2), MPI_DOUBLE, 0, 0, newComm, MPI_STATUS_IGNORE);
    MPI_Recv(&(B[1][1]), (N-2)*(N-2), MPI_DOUBLE, 0, 0, newComm, MPI_STATUS_IGNORE);
    
    printf("Recive subarray size %d \n", (N-2)*(N-2));
    

    value_t ghost_left[N];
    value_t ghost_right[N];
    value_t ghost_up[N];
    value_t ghost_down[N];


	int left_rank;
	int right_rank;
	int up_rank;
	int down_rank;
	MPI_Cart_shift(newComm, 0, 1, &left_rank, &right_rank);
	MPI_Cart_shift(newComm, 1, 1, &up_rank, &down_rank);
    
    //for each time step
    time_t start = time(NULL);
    for (int t = 0; t < T; t++) {
		
		MPI_Send(&A[1], N, MPI_DOUBLE, up_rank, 0, newComm);
		
		MPI_Send(&A[N-2], N, MPI_DOUBLE, down_rank, 0, newComm);
		
		value_t *tempArray = getColumn(1,A,N);
		MPI_Send(&tempArray, N, MPI_DOUBLE, left_rank, 0, newComm);
		
		tempArray = getColumn(N-2,A,N);
		MPI_Send(&tempArray, N, MPI_DOUBLE, right_rank, 0, newComm);
		
		free(tempArray);
		
		MPI_Recv(&(ghost_up[0]), N, MPI_DOUBLE, up_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(&(ghost_down[0]), N, MPI_DOUBLE, down_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(&(ghost_left[0]), N, MPI_DOUBLE, left_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(&(ghost_right[0]), N, MPI_DOUBLE, right_rank, 0, newComm, MPI_STATUS_IGNORE);
		
		A[0] = ghost_up;
		A[N-1] = ghost_down;
		B[0] = ghost_up;
		B[N-1] = ghost_down;
		insertColumn(ghost_left, 0, A, N);
		insertColumn(ghost_left, 0, B, N);
		insertColumn(ghost_right, N, A, N);
		insertColumn(ghost_right, N, B, N);

		
        //we propagate the temparature
        for (long long i = 0; i < N; i++) {
            for (long long j = 0; j < N; j++) {
                if (i == Y && j == X) {
                    B[i][j] = A[i][j];
                    continue;
                }

                //get temperatur at current position
                value_t tc = A[i][j];

                //get temperatur of adjacent cells
                value_t t_above = (i != 0) ? A[i - 1][j] : tc;
                value_t t_left = (j != 0) ? A[i][j - 1] : tc;
                value_t t_right = (j != N - 1) ? A[i][j + 1] : tc;
                value_t t_below = (i != N - 1) ? A[i + 1][j] : tc;

                B[i][j] = tc + 0.2 * (
                        t_above + t_left + 
                        t_right + t_below + ( -4 * tc )
                    );
            }
        }

        Vector H = A;
        A = B;
        B = H;

		if(rank == 0) {
			if (!(t % 1000)) {
				printf("Current timestamp t=%d\n", t);
			}
		}
    }
    time_t stop = time(NULL);

    releaseVector(B, N);

    printf("Verification from rank %d: %s\n", rank, (verify(A, N)) ? "OK" : "FAILED");
    
    
    //release the vector again
    releaseVector(A, N);

	if(rank == 0) {

		printf("The Program took %ld to execute\n", stop - start);
	}

    
    MPI_Finalize();
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
        //free(m[i]);
    }

    //free(m);
}

void fill_vector(Vector m, int size, int x, int y) {
    if (x > size || y > size) {
        return;
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            m[i][j] = 273;
        }
    }

    m[y][x] = 273 + 60;

}

int verify(Vector m, int size) {
    int success = 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            value_t temp = m[i][j];
            if (273 <= temp && temp <= 273 + 60) {
                continue;
            }
            printf("%f --> x: %d y: %d\n", temp, j, i);
            success = 0;
            break;
        }
    }
    return success;
}

value_t *getColumn(int pos, Vector m, int size) {
	
	value_t temp[size];
	
	for(int i = 0; i < size; i++) {
		temp[i] = m[i][pos];
	}
	
	return temp;
	
}


void insertColumn(value_t *c, int pos, Vector m, int size) {
	
	for(int i = size-2; i >= 0; i--) {
		for(int j = size-1; j >= 0; j--) {
			m[i][j+1] = m[i][j];
			if(j == pos) {
				m[i][j] = c[j];
				break;
			}
		}
	}
		
	
}

