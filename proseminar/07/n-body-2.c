#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

#define G 1
#define SEED 55754186

typedef struct {
	int x;
	int y;
	double velocity_x;
	double velocity_y;
	int mass;
	int identifier;
} Particle;

Particle* create_particles(int number_of_particles);
void init_particles(Particle* particles, int number_of_particles);
void release_particles(Particle *particles);
void print_particles(Particle *particles, int number_of_particles);

double calcForce(Particle p1, Particle p2);
double distance(Particle p1, Particle p2);
double calcVelocity(double force, Particle particle);
Particle updatePostion(Particle particle, int size);
double sum_up_force(Particle *particles, int number_of_particles, int index);
Particle* calculate_new_timestamp(Particle *particles, Particle *allParticles, int number_of_particles, int number_of_all_particles);
Particle *calculate_new_timestamp_para(Particle *particles, Particle *allParticles, int number_of_particles, int number_of_all_particles);

void swap(int *a, int *b);
void randomize(int arr[], int n);

int main(int argc, char **argv) {
    if (argc < 3) {
		printf("First parameter is the number of particles for each rank, second parameter is the number of timestamps\n");
		return EXIT_FAILURE;
	}

	//int spaceSize = atoi(argv[1]);
	int numberParticles = atoi(argv[1]);
	int timestamps = atoi(argv[2]);

    int numProcs = 8;

    Particle* particles = create_particles(numberParticles * numProcs);
    
	init_particles(particles, numberParticles * numProcs);

    struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);

    omp_set_num_threads(numProcs);
    

	for (int i = 0; i < timestamps; i++) {
        Particle *temp = calculate_new_timestamp(particles, particles, numberParticles * numProcs, numberParticles * numProcs);
        release_particles(particles);
        particles = temp;
	}

    
    gettimeofday(&tv2, NULL);

    printf ("%d;%f;\n",numberParticles * numProcs,
            (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
            (double) (tv2.tv_sec - tv1.tv_sec));
    

    release_particles(particles);

    return EXIT_SUCCESS;
}

Particle* create_particles(int number_of_particles) {
    Particle* particles = (Particle*)malloc(sizeof(Particle) * number_of_particles);
    if (!particles) {
        perror("Error allocating memory");
        return NULL;
    }

    return particles;
}

void init_particles(Particle* particles, int number_of_particles) {
	int m;

    int steps = (int)(INT_MAX / number_of_particles);

    int *random_x_values = malloc(sizeof(int)* number_of_particles);
    if (!random_x_values) {
        perror("Error allocating memory");
        return;
    }
    int *random_y_values = malloc(sizeof(int)* number_of_particles);
    if (!random_y_values) {
        perror("Error allocating memory");
        return;
    }

    srand(time(NULL));
    for (int i = 0; i < number_of_particles; i++) {
        random_x_values[i] = (rand() % steps) + (i * steps);
        random_y_values[i] = (rand() % steps) + (i * steps);
    }

    randomize(random_x_values, number_of_particles);
    randomize(random_y_values, number_of_particles);

    for (int i = 0; i < number_of_particles; i++) {
        Particle particle;
		m = 1;

		particle.velocity_x = rand() % 100;//so v doesn't get too high
		particle.velocity_y = rand() % 100;
		particle.mass = m;

		particle.x = random_x_values[i];
		particle.y = random_y_values[i];

		particle.identifier = i + 1;
        particles[i] = particle;
    }

    free(random_x_values);
    free(random_y_values);
}

void print_particles(Particle *particles, int number_of_particles) {
    for (int i = 0; i < number_of_particles; i++) {
        printf("x: %d\t--\ty: %d \tm: %d\n", particles[i].x, particles[i].y, particles[i].mass);
    }

    printf("\n");
}

void release_particles(Particle* particles) {
    free(particles);
}


double calcForce(Particle p1, Particle p2) {

	double r = distance(p1, p2);
	return G * (p1.mass * p2.mass) / r;

}

double distance(Particle p1, Particle p2) {
	return pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2);
}

double calcVelocity(double force, Particle particle) {
	return (force / particle.mass);

}

Particle updatePostion(Particle particle, int size) {

	long x = particle.x;
	long y = particle.y;
	double v_x = particle.velocity_x;
	double v_y = particle.velocity_x;

	x = x + v_x;
	x %= size;

	y = y + v_y;
	y %= size;

	//printf("vx%f vy%f\n", particle.velocity_x, particle.velocity_y);

	//printf("oldx%d newx%d\n", particle.x, x);
	//printf("oldy%d newy%d\n", particle.y, y);

	particle.x = x;
	particle.y = y;

	return particle;

}

double sum_up_force(Particle *particles, int number_of_particles, int index) {
    double force = 0.0;
    for (int i = 0; i < number_of_particles; i++) {
        if (i != index) {
            force += calcForce(particles[index], particles[i]);
        }
    }
    return force;
}

Particle *calculate_new_timestamp(Particle *particles, Particle *allParticles, int number_of_particles, int number_of_all_particles) {
    Particle *temp = create_particles(number_of_particles);

    #pragma omp for nowait
    for (int i = 0; i < number_of_particles; i++) {
    	temp[i] = particles[i];
    	for(int j = 0; j < number_of_all_particles; j++) {
    		if(temp[i].identifier != allParticles[j].identifier) {
    			//printf("m1%d m2%d\n", temp[i].mass, allParticles[i].mass);

				double force = calcForce(temp[i], allParticles[j]);
				//https://gamedev.stackexchange.com/questions/48119/how-do-i-calculate-how-an-object-will-move-from-one-point-to-another
				double angle = atan2(temp[i].y - allParticles[j].y,  temp[i].x - allParticles[j].x);
				temp[i].velocity_x += calcVelocity(force, allParticles[j]) * cos(angle);
				temp[i].velocity_y += calcVelocity(force, allParticles[j]) * sin(angle);

    		}
    	}
    	 temp[i] = updatePostion(temp[i], INT_MAX);
    }
    #pragma omp barrier
    
    return temp;
}

Particle *calculate_new_timestamp_para(Particle *particles, Particle *allParticles, int number_of_particles, int number_of_all_particles) {
    Particle *temp = create_particles(number_of_particles);

    //#pragma omp for 
    for (int i = 0; i < number_of_particles; i++) {
        temp[i] = particles[i]; 
        for(int j = 0; j < number_of_all_particles; j++) {
            if(temp[i].identifier != allParticles[j].identifier) {
                //printf("m1%d m2%d\n", temp[i].mass, allParticles[i].mass);

                double force = calcForce(temp[i], allParticles[j]);
                //https://gamedev.stackexchange.com/questions/48119/how-do-i-calculate-how-an-object-will-move-from-one-point-to-another
                double angle = atan2(temp[i].y - allParticles[j].y,  temp[i].x - allParticles[j].x);
                temp[i].velocity_x += calcVelocity(force, allParticles[j]) * cos(angle);
                temp[i].velocity_y += calcVelocity(force, allParticles[j]) * sin(angle);

            }
        }
        
        temp[i] = updatePostion(temp[i], INT_MAX);
    }
    
    
    return temp;
}

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void randomize(int arr[], int n) {
    srand(SEED);
    int i;
    for(i = n-1; i > 0; i--) {
        int j = rand() % (i+1);
        swap(&arr[i], &arr[j]);
    }
}
