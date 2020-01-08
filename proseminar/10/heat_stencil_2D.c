#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <omp.h>

typedef double value_t;

// -- vector utilities --

typedef value_t **Vector;

Vector createVector(int N, int M);

void releaseVector(Vector m, int size);

void fill_vector(Vector m, int size, int x, int y);

int verify(Vector m, int size);

void printVector(Vector m, int size);

value_t *getColumn(int pos, Vector m, int size);

void insertColumn(value_t *c, int pos, Vector m, int size);

// -- simulation code ---

int main(int argc, char **argv) {
    //width of stencil field
    int N_big = 10; 

    //number of timestamps
    int T = 1000;
    if (argc > 1) {
        N_big = atoi(argv[1]);
    }
    
    /* initialize mpi*/
	int rank;
	int numProcs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    //start calculation----------------------------------------------------------------------------------
	//+2 because left and right ghost cell
    int temp_val = ceil((double)N_big/numProcs);
	int N = rank == numProcs - 1 && N_big % temp_val != 0 ? N_big % temp_val : temp_val;

	//the position of the heat-source
	int X = N_big / 4;
	int Y = N_big / 5;


    //recieve send data
    int vec_size = N;
    Vector A = createVector(vec_size + 2, N_big);
    Vector B = createVector(vec_size + 2, N_big);

    //create the stencil field
	if(rank == 0) {
        //create a buffer
		Vector A_big = createVector(N_big, N_big);
		fill_vector(A_big, N_big, X, Y);
        
        //copy for rank 0
        for (int i = 1; i < vec_size + 1; i++) {
            for (int j = 0; j < N_big; j++) {
                A[i][j] = A_big[i - 1][j];
                
            }  
        }
        
		for(int i = 1; i < numProcs; i++) {
            int vec_size = i == numProcs - 1 && N_big % temp_val != 0 ? N_big % temp_val : temp_val; 
            MPI_Send(&(A_big[i * N][0]), N_big * vec_size, MPI_DOUBLE, i, 42, MPI_COMM_WORLD);
            	    
		}
			
		releaseVector(A_big, N_big);
	}
	
    if (rank != 0) {
	    MPI_Recv(&(A[1][0]), (N_big)*(vec_size), MPI_DOUBLE, 0, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
	
    
    //for each time step
    MPI_Request send_above_request;
    MPI_Request send_below_request;

    
    time_t start = time(NULL);
    for (int t = 0; t < T; t++) {
		
		
        //send line below
        if (rank != numProcs - 1) {
            MPI_Isend(&(A[vec_size][0]), N_big, MPI_DOUBLE, rank + 1, 41, MPI_COMM_WORLD, &send_below_request);
        }

        //get line above
        if (rank != 0) {
            MPI_Irecv(&(A[0][0]), N_big, MPI_DOUBLE, rank - 1, 41, MPI_COMM_WORLD, &send_below_request);
        }
        
        //send line above
        if (rank != 0) {
            MPI_Isend(&(A[1][0]), N_big, MPI_DOUBLE, rank - 1, 40, MPI_COMM_WORLD, &send_above_request);
        }


        //get line below
        if (rank != numProcs - 1) {
            MPI_Irecv(&(A[vec_size + 1][0]), N_big, MPI_DOUBLE, rank + 1, 40, MPI_COMM_WORLD, &send_above_request);
        }

        MPI_Wait(&send_below_request, MPI_STATUS_IGNORE);
        MPI_Wait(&send_above_request, MPI_STATUS_IGNORE);

        
        #pragma omp parallel for schedule(guided) shared(A, B) firstprivate(N, Y, X, vec_size, rank, N_big, numProcs) default(none) 
        for (long long i = 1; i < vec_size + 1; i++) {
            for (long long j = 0; j < N_big; j++) {
                if (i - 1 == Y && j == X) {
                    B[i][j] = A[i][j];
                    continue;
                }

                //get temperatur at current position
                value_t tc = A[i][j];

                //get temperatur of adjacent cells
                value_t t_above = (rank != 0 || i != 1) ? A[i - 1][j] : tc;
                value_t t_left = (j != 0) ? A[i][j - 1] : tc;
                value_t t_right = (j != N_big - 1) ? A[i][j + 1] : tc;
                value_t t_below = (rank != numProcs - 1 || i != vec_size) ? A[i + 1][j] : tc;

                B[i][j] = tc + 0.2 * (
                        t_above + t_left + 
                        t_right + t_below + ( -4 * tc )
                    );
            }
        }

        Vector H = A;
        A = B;
        B = H;
    }

    //gather everything
    Vector A_big = createVector(N_big, N_big);
    
    MPI_Request *gath = malloc((numProcs - 1) * sizeof(MPI_Request));
    if (!gath) {
        perror("Could not allocate memory");
        releaseVector(A_big, N_big);
        releaseVector(B, vec_size + 2);
        releaseVector(A, vec_size + 2);
        MPI_Finalize();

        return EXIT_SUCCESS;
    }

    if (rank != 0) {
        MPI_Isend(&(A[1][0]), N_big * vec_size, MPI_DOUBLE, 0, 30, MPI_COMM_WORLD, &(gath[rank - 1]));
    } else {
        for (int i = 1; i < vec_size + 1;i++) {
            for (int j = 0; j < N_big; j++) {
                A_big[i-1][j] = A[i][j];
            }
        }

        for (int i = 1; i < numProcs; i++) {
            MPI_Irecv(&(A_big[i * vec_size][0]), N_big * vec_size, MPI_DOUBLE, i, 30, MPI_COMM_WORLD,&(gath[i - 1]));
            MPI_Wait(&(gath[i - 1]), MPI_STATUS_IGNORE);
        }
    }

    if (rank == 0) {
        printf("Result: \n");
       for (int i = 0;i < N_big; i++) {
            for (int j = 0; j < N_big; j++) {
                printf("%f (%d-%d: %d)\t", A_big[i][j],i,j, rank);
            }
            printf("\n");
        }
    }


    free(gath);
    releaseVector(A_big, N_big);
    releaseVector(B, vec_size + 2);
    releaseVector(A, vec_size + 2);
    MPI_Finalize();

    return EXIT_SUCCESS;
}

Vector createVector(int N, int M) {
    value_t *data = calloc(N * M, sizeof(value_t));
    Vector arr = calloc(N, sizeof(value_t *));

    if (data == NULL || arr == NULL) {
        perror("Could not allocate memory");
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        arr[i] = &(data[i*M]);
    }

    return arr;
}

void releaseVector(Vector m, int size) {
    free(m[0]);
    free(m);
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
	
	value_t *temp = (value_t*)malloc(size * sizeof(value_t));
    if (!temp) {
        perror("Could not allocate memory");
        return NULL;
    }
	
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

