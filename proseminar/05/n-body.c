#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 1
#define SEED 55754186

typedef struct {
	int x;
	int y;
	double velocity;
	int mass;
	int identifier;
} Particle;

typedef Particle **Space;

Space createSpace(int size);
void initSpace(Space space, int size, int number);
void printSpace(Space space, int size, char *title);
double calcForce(Particle p1, Particle p2);
double distance(Particle p1, Particle p2);
double calcVelocity(double force, Particle particle);
Space calcNewSpace(Space space, int size);
Particle updatePostion(Particle particle, int size);
void releaseSpace(Space space);

int main(int argc, char **argv) {
	if (argc < 4) {
		printf(
				"First parameter is size of the quadratic space, second parameter is number of particles and the third is the number of timestamps\n");
		return EXIT_FAILURE;
	}

	int spaceSize = atoi(argv[1]);
	int numberParticles = atoi(argv[2]);
	int timestamps = atoi(argv[3]);

	if (numberParticles > spaceSize * spaceSize) {
		printf("Number of particles must be smaller than space size (N*N)\n");
		return EXIT_FAILURE;
	}

	Space space = createSpace(spaceSize);

	initSpace(space, spaceSize, numberParticles);

	printSpace(space, spaceSize, "Space");

	for (int i = 0; i < timestamps; i++) {
		space = calcNewSpace(space, spaceSize);
		sleep(1);

		char *timestamp_string = malloc(sizeof(char) * 20);
		if (!timestamp_string) {
			perror("Could not alloc memory");
		}

		sprintf(timestamp_string, "New Space %d", i + 1);
		printSpace(space, spaceSize, timestamp_string);
		free(timestamp_string);
	}

	releaseSpace(space);

	return EXIT_SUCCESS;
}

/*
 * creates a size*size array
 */
Space createSpace(int size) {
	Particle *data = (Particle*) calloc(size * size, sizeof(Particle));
	Space array = (Particle**) calloc(size, sizeof(Particle*));
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

		particle.identifier = i + 1;
		space[x][y] = particle;

	}

}

void printSpace(Space space, int size, char *title) {

	printf("\033[2J");
	printf("%s\n", title);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (space[i][j].mass == 0) {
				printf("-\t");
			} else {
				printf("%.2f (%d)\t", space[i][j].velocity, space[i][j].identifier);
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
	return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));

}

double calcVelocity(double force, Particle particle) {
	return particle.velocity + force / particle.mass;

}

/*
 * calculates the velocity of all particles and updates the position
 */
Space calcNewSpace(Space space, int size) {

	Space newSpace = createSpace(size);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			//empty space
			if (space[i][j].mass == 0)
				continue;

			Particle temp = space[i][j];
			double v = 0;

			for (int k = 0; k < size; k++) {
				for (int l = 0; l < size; l++) {
					//empty space or same
					if (space[k][l].mass == 0 || (k == i && l == j))
						continue;
					double force = calcForce(temp, space[k][l]);
					v += calcVelocity(force, temp);
				}
			}

			temp.velocity = v;
			//printf("velocity%f\n", temp.velocity);
			Particle newParticle = updatePostion(temp, size);

			//colliding
			if (newSpace[newParticle.x][newParticle.y].mass > 0)
				newSpace[newParticle.x][newParticle.y].velocity += newParticle.velocity;
			else
				newSpace[newParticle.x][newParticle.y] = newParticle;

		}
	}

	releaseSpace(space);
	return newSpace;

}

Particle updatePostion(Particle particle, int size) {

	int x = particle.x;
	int y = particle.y;
	double v = particle.velocity;

	x = x + v;
	x %= size;

	y = y + v;
	y %= size;

	//printf("oldx%d newx%d\n", particle.x, x);
	//printf("oldy%d newy%d\n", particle.y, y);

	particle.x = x;
	particle.y = y;

	return particle;

}

void releaseSpace(Space space) {
	free(space[0]);
	free(space);
}

