#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

// -- vector utilities --

typedef value_t ***Vector;

Vector createVector(int N);

void releaseVector(Vector m, int size);

void fill_vector(Vector m, int size, int x, int y, int z);

int verify(Vector m, int size);

value_t *getColumn(int pos, Vector m, int size);

void columnShiftRight(int start, Vector m, int size);

void insertColumn(value_t *c, int pos, Vector m, int size);

// -- simulation code ---

int main(int argc, char **argv) {
    int N_big = 10;
    if (argc > 1) {
        N = atoi(argv[1]);
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

    int T = N * 500;
    
    //+4 because left and right behind before ghost cell
	int N = (N_big/numProcs) + 4;
	
	Vector A = createVector(N);
	Vector B = createVector(N);

	//fill the array
	int X = N_big / 4;
	int Y = N_big / 5;
	int Z = N_big / 2;
	fill_vector(A, N, X, Y, Z);
    
    if(rank == 0) {
		printf("Computing heat-distribution for romm size %dx%d for %d timestamps\n", N_big, N_big, T);

		//create a buffer
		Vector A_big = createVector(N_big);
		
		fill_vector(A_big, N_big, X, Y, Z);

		//create a second buffer for the computation
		Vector B_big = createVector(N_big);
		
		for(int i = 0; i < numProcs; i++) {
				
			int array_size[] = {N_big,N_big, N_big};
			int array_subsize[] = {N-4,N-4,N-4};
			int array_start[] = {0,0,0};
			
			MPI_Datatype subArray;

			MPI_Type_create_subarray(3, array_size, array_subsize, array_start, MPI_ORDER_C, MPI_INT, &subArray);
			MPI_Type_commit(&subArray);
			
			MPI_Send(&(A_big[0][0][0]), 1, subArray, i, 0, MPI_COMM_WORLD);
			MPI_Send(&(B_big[0][0][0]), 1, subArray, i, 0, MPI_COMM_WORLD);
			
			MPI_Type_free(&subArray);
				
		}
	}
	
	//A[1][1][1] - first cell left out for ghost cells
	MPI_Recv(&(A[1][1][1]), (N-4)*(N-4)*(N-4), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(B[1][1][1]), (N-4)*(N-4)*(N-4), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    printf("Recive subarray size %d \n", (N-4)*(N-4));
    
	

    value_t **ghost_left;
    value_t **ghost_right;
    value_t **ghost_up;
    value_t **ghost_down;
    value_t **ghost_behind;
    value_t **ghost_before;
    
    MPI_Comm newComm;
    
    int dims[] = {N,N,N};
    int periods[] = {1,1,1};
    int reorder = 1;
    
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &newComm);

	int left_rank;
	int right_rank;
	int up_rank;
	int down_rank;
	int behind_rank;
	int before_rank
	MPI_Cart_shift(newComm, 0, 1, &left_rank, &right_rank);
	MPI_Cart_shift(newComm, 1, 1, &up_rank, &down_rank);
	MPI_Cart_shift(newComm, 2, 1, &before_rank, &behind_rank);
    
    //for each time step
    for (int t = 0; t < T; t++) {
		
		MPI_Send(&A[1], N, MPI_DOUBLE, up_rank, 0, newComm);
		MPI_Send(&A[N-2], N, MPI_DOUBLE, down_rank, 0, newComm);
		value_t **tempArray = getCell(1,1,A,N);
		MPI_Send(&tempArray, N, MPI_DOUBLE, left_rank, 0, newComm);
		tempArray = getCell(N-1,1,A,N);
		MPI_Send(&tempArray, N, MPI_DOUBLE, right_rank, 0, newComm);
		tempArray = getCell(1,N-1,A,N);
		MPI_Send(&tempArray, N, MPI_DOUBLE, behind_rank, 0, newComm);
		tempArray = getCell(N-1,N-1,A,N);
		MPI_Send(&tempArray, N, MPI_DOUBLE, before_rank, 0, newComm);
		
		MPI_Recv(ghost_up, 1, MPI_DOUBLE, up_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_down, 1, MPI_DOUBLE, down_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_left, 1, MPI_DOUBLE, left_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_right, 1, MPI_DOUBLE, right_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_behind, 1, MPI_DOUBLE, behind_rank, 0, newComm, MPI_STATUS_IGNORE);
		MPI_Recv(ghost_before, 1, MPI_DOUBLE, before_rank, 0, newComm, MPI_STATUS_IGNORE);
		
		A[0] = ghost_up;
		A[N-1] = ghost_down;
		B[0] = ghost_up;
		B[N-1] = ghost_down;
		insertZell(ghost_left, 0,0, A, N);
		insertZell(ghost_left, 0,0, B, N);
		insertZell(ghost_right, N,0, A, N);
		insertZell(ghost_right, N,0, B, N);
		insertZell(ghost_behind, 0,N, A, N);
		insertZell(ghost_behind, 0,N, B, N);
		insertZell(ghost_before, N,N, A, N);
		insertZell(ghost_before, N,N, B, N);
		
		
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


                    B[i][j][k] = tc + 0.2 * (
                              t_above + t_left + 
                        t_right + t_below + + t_behind + t_before + ( -6 * tc )
                        );
                }
            }
        }

        Vector H = A;
        A = B;
        B = H;

        if (!(t % 1000)) {
            printf("Current timestamp t=%d\n", t);
        }
    }

    releaseVector(B, N);

    printf("Verification: %s\n", (verify(A, N)) ? "OK" : "FAILED");

    //release the vector again
    releaseVector(A, N);
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

void releaseVector(Vector m, int size) {

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            free(m[i][j]);
        }
        free(m[i]);
    }

    free(m);
}

void fill_vector(Vector m, int size, int x, int y, int z) {
    if (x > size || y > size) {
        return;
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                m[i][j][k] = 273;
            }
        }
    }

    m[y][x][z] = 273 + 60;

}

int verify(Vector m, int size) {
    int success = 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                value_t temp = m[i][j][k];
                if (273 <= temp && temp <= 273 + 60) {
                    continue;
                }
                printf("%f --> x: %d y: %d z: %d\n", temp, j, i, k);
                success = 0;
                break;
            }
        }
    }
    return success;
}


//TODO
value_t *getZell(int y, int z, Vector m, int size) {
	
	value_t temp[size][size];
	
	for(int i = 0; i < size; i++) {
		temp[i] = m[i][y];
	}
	
	return temp;
	
}

void zellShiftRight(int start, Vector m, int size) {
		
	for(int i = size; i <= start; i++) {
		m[i] = getColumn(i+1,m,size);
	}
	
}


void insertZell(value_t *c, int y, int z, Vector m, int size) {
	
	for (int i = 0; i < size; i++) {
		if(i == pos) {
			columnShiftRight(i,m,size);
			m[i] = c;
		}
	}
	
	
}
