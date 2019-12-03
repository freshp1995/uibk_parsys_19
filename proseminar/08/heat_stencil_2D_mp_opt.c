#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>

typedef double value_t;

// -- vector utilities --

typedef value_t **Vector;

Vector createVector(int N);

void releaseVector(Vector m, int size);

void fill_vector(Vector m, int size, int x, int y);

int verify(Vector m, int size);

// -- simulation code ---
long timediff(clock_t t1, clock_t t2) {
    long elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000;
    return elapsed;
}

int main(int argc, char **argv) {
    int N = 50;
    if (argc > 1) {
        N = atoi(argv[1]);
    }

    int T = 5000;
    //printf("Computing heat-distribution for romm size %dx%d for %d timestamps\n", N, N, T);

    //create a buffer
    Vector A = createVector(N);

    //fill the array
    int X = N / 4;
    int Y = N / 5;
    fill_vector(A, N, X, Y);


    //create a second buffer for the computation
    Vector B = createVector(N);
    
    struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);
    //for each time step
    for (int t = 0; t < T; t++) {
        //we propagate the temparature
        #pragma omp parallel for
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

        if (!(t % 1000)) {
            //printf("Current timestamp t=%d\n", t);
        }
    }

    releaseVector(B, N);

    gettimeofday(&tv2, NULL);

    printf ("%f;\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));

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
    }

    return y;
}

void releaseVector(Vector m, int size) {
    for (int i = 0; i < size; i++) {
        free(m[i]);
    }

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
            //printf("%f --> x: %d y: %d\n", temp, j, i);
            success = 0;
            break;
        }
    }
    return success;
}
