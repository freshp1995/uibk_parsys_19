#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

// -- vector utilities --

typedef value_t ***Vector;

Vector createVector(int N);

void releaseVector(Vector m, int size);

void fill_vector(Vector m, int size, int x, int y, int z);

int verify(Vector m, int size);

// -- simulation code ---

int main(int argc, char **argv) {
    int N = 10;
    if (argc > 1) {
        N = atoi(argv[1]);
    }

    int T = N * 500;
    printf("Computing heat-distribution for romm size %dx%d for %d timestamps\n", N, N, T);

    //create a buffer
    Vector A = createVector(N);

    //fill the array
    int X = N / 4;
    int Y = N / 5;
    int Z = N / 2;
    fill_vector(A, N, X, Y, Z);


    //create a second buffer for the computation
    Vector B = createVector(N);
    
    //for each time step
    for (int t = 0; t < T; t++) {
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
