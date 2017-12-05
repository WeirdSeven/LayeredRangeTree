#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

double randfrom(double min, double max) {
	return min + (((double) rand()) / RAND_MAX * (max - min));
}

int main() {
	long long unsigned num_points;

	//printf("Enter number of points: ");
	scanf("%llu", &num_points);
	printf("%llu\n", num_points);

	if ((num_points & (num_points - 1)) != 0) {
		printf("Please enter a power of 2.\n");
		exit(EXIT_FAILURE);
	}

	double square_root = sqrt(num_points) / 2;

	srand(time(NULL));

	for(long long unsigned i = 0; i < num_points; i++) {
		double x = randfrom(-square_root, square_root);
		double y = randfrom(-square_root, square_root);
		printf("%f %f\n", x, y);
	}

	return 0;
}