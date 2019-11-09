#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 1
#define SEED 55754186

typedef struct {
	int x;
	int y;
	double velocity;
	int mass;
} Particle;

typedef Particle **Space;

Space createSpace(int size);
void initSpace(Space space, int size, int number);
void printSpace(Space space, int size);
double calcForce(Particle p1, Particle p2);
double distance(Particle p1, Particle p2);
double calcVelocity(double force, Particle particle);
void calcNewSpace(Space space, int size);
Particle updatePostion(Particle particle, int size);

int main(int argc, char **argv) {
	if (argc < 3) {
			printf("First parameter is size of the quadratic space, second parameter is number of particles\n");
			return EXIT_FAILURE;
		}


	int spaceSize = atoi(argv[1]);
	int numberParticles = atoi(argv[2]);

	if(numberParticles > spaceSize*spaceSize) {
		printf("Number of particles must be smaller than space size (N*N)\n");
		return EXIT_FAILURE;
	}

	Space space = createSpace(spaceSize);

	initSpace(space, spaceSize, numberParticles);

	printSpace(space, spaceSize);

	printf("\n\n");
	calcNewSpace(space, spaceSize);

	printSpace(space, spaceSize);

}

/*
 * creates a size*size array
 */
Space createSpace(int size) {
	Particle *data = (Particle*) malloc(size * size * sizeof(Particle));
	Space array = (Particle**) malloc(size * sizeof(Particle*));
	for (int i = 0; i < size; i++)
		array[i] = &(data[size * i]);

	return array;

}

/*
 * initialize the space with number particles
 */
void initSpace(Space space, int size, int number) {
	srand(SEED);

	double v;
	int x, y, m;

	for (int i = 0; i < number; i++) {

		Particle particle;
		v = rand() % size; //so v doesn't get too high
		m = 1;

		particle.velocity = v;
		particle.mass = m;

		// generate x and y until we find empty space
		do {
			x = rand() % size;
			y = rand() % size;

			particle.x = x;
			particle.y = y;

		//if no particle there, then mass is 0
		} while (space[x][y].mass > 0);

		space[x][y] = particle;

	}

}

void printSpace(Space space, int size) {

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (space[i][j].mass == 0) {
				printf("-\t");
			} else {
				printf("%f\t", space[i][j].velocity);
			}
		}
		printf("\n");
	}

}

double calcForce(Particle p1, Particle p2) {

	double r = distance(p1, p2);
	return G * (p1.mass * p2.mass) / pow(r, 2);

}

double distance(Particle p1, Particle p2) {
	return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2) * 1.0);

}

double calcVelocity(double force, Particle particle) {
	return particle.velocity + force / particle.mass;

}

/*
 * calculates the velocity of all particles and updates the postion
 */
void calcNewSpace(Space space, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (space[i][j].mass > 0) {
				Particle temp = space[i][j];

				for (int k = 0; k < size; k++) {
					for (int l = 0; l < size; l++) {
						double force = calcForce(temp, space[k][l]);
						//printf("force%f\n", force);
						temp.velocity = calcVelocity(force, temp);
					}
				}
				temp = updatePostion(temp, size);

			}
		}
	}
}

Particle updatePostion(Particle particle, int size) {

	int x = particle.x;
	int y = particle.y;
	double v = particle.velocity;

	x = x + v;
	x %= size;

	y = y + v;
	y %= size;

	printf("oldx%d newx%d\n", particle.x, x);

	particle.x = x;
	particle.y = y;


	return particle;

}
