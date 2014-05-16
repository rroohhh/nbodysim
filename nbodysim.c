/*
 * nbodysim.c
 *
 * An fast n-body simulation
 *
 * Created on: 09.05.2014
 * Author: robin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

//Global Vars
double t;
double gc = 6.67384e-11;
int size;
double time;
double endtime;
char output[100];
FILE *file;

//Typedefs:

//vector
typedef struct {
	double x;
	double y;
	double z;
} vector;

//matrix
typedef struct {
	vector old;
	vector now;
	vector new;
} matrix;

//body
typedef struct {
	vector v;
	matrix r;
	double mass;
	vector f;
} body;

body *univers = NULL;

void sum_grav(int body) {
	int i;
	double tmp;
	double mgrav;

	for (i = body + 1; i < size; i++) {
		double dx = (univers[i].r.now.x - univers[body].r.now.x);
		double dy = (univers[i].r.now.y - univers[body].r.now.y);
		double dz = (univers[i].r.now.z - univers[body].r.now.z);
		double dist = sqrt(dx * dx + dy * dy + dz * dz);
		double diff = dist * dist * dist;
		/*
		 printf(
		 "dist(%i)(x)(y)(z)=%15.10lf,%15.10lf,%15.10lf,gesamtdist: %15.10lf,diff(dist^3): %15.10lf\n",
		 i, dx, dy, dz, dist, diff);
		 */
		mgrav = gc * univers[i].mass * univers[body].mass;
		//printf("mgrav: %15.10lf", mgrav);
		//printf("Masse: %15.10lf, %15.10lf\n", univers[body].mass,	univers[i].mass);
		tmp = mgrav * (dx / diff);
		univers[body].f.x += tmp;
		univers[i].f.x += -1 * tmp;
		tmp = mgrav * (dy / diff);
		univers[body].f.y += tmp;
		univers[i].f.y += -1 * tmp;
		tmp = mgrav * (dz / diff);
		univers[body].f.z += tmp;
		univers[i].f.z += -1 * tmp;
	}
	//printf("%i. Körper %15.10lf,%15.10lf,%15.10lf\n", body, univers[body].f.x,
	//		univers[body].f.y, univers[body].f.z);
}

void leap(int body) {
	/*
	 printf("pos_old(%i)(x)(y)(z)=%15.10lf,%15.10lf,%15.10lf\n", body,
	 univers[body].r.old.x, univers[body].r.old.y,
	 univers[body].r.old.z);
	 printf("pos_now(%i)(x)(y)(z)=%15.10lf,%15.10lf,%15.10lf\n", body,
	 univers[body].r.now.x, univers[body].r.now.y,
	 univers[body].r.now.z);
	 */
	univers[body].r.new.x = 2 * univers[body].r.now.x - univers[body].r.old.x
			+ (univers[body].f.x / univers[body].mass) * t * t;
	univers[body].r.new.y = 2 * univers[body].r.now.y - univers[body].r.old.y
			+ (univers[body].f.y / univers[body].mass) * t * t;
	univers[body].r.new.z = 2 * univers[body].r.now.z - univers[body].r.old.z
			+ (univers[body].f.z / univers[body].mass) * t * t;

	univers[body].r.old.x = univers[body].r.now.x;
	univers[body].r.old.y = univers[body].r.now.y;
	univers[body].r.old.z = univers[body].r.now.z;
	univers[body].r.now.x = univers[body].r.new.x;
	univers[body].r.now.y = univers[body].r.new.y;
	univers[body].r.now.z = univers[body].r.new.z;

}

void setup_next() {
	int i;
	for (i = 0; i < size; i++) {
		univers[i].f.x = 0;
		univers[i].f.y = 0;
		univers[i].f.z = 0;
	}
}

void setup_univers() {
	int i;

	for (i = 0; i < size; i++) {
		univers[i].f.x = 0;
		univers[i].f.y = 0;
		univers[i].f.z = 0;
	}

	for (i = 0; i < size; i++) {
		sum_grav(i);
	}

	for (i = 0; i < size; i++) {
		univers[i].r.new.x = univers[i].r.now.x + univers[i].v.x * t
				+ 0.5 * (univers[i].f.x / univers[i].mass) * t * t;
		univers[i].r.new.y = univers[i].r.now.y + univers[i].v.y * t
				+ 0.5 * (univers[i].f.y / univers[i].mass) * t * t;
		univers[i].r.new.z = univers[i].r.now.z + univers[i].v.z * t
				+ 0.5 * (univers[i].f.z / univers[i].mass) * t * t;
		univers[i].r.old.x = univers[i].r.now.x;
		univers[i].r.old.y = univers[i].r.now.y;
		univers[i].r.old.z = univers[i].r.now.z;
		univers[i].r.now.x = univers[i].r.new.x;
		univers[i].r.now.y = univers[i].r.new.y;
		univers[i].r.now.z = univers[i].r.new.z;
	}
}

void calc_next() {
	int i;

	setup_next();

	for (i = 0; i < size; i++) {
		sum_grav(i);
	}

	for (i = 0; i < size; i++) {
		leap(i);
	}
}

void user_prompt() {
	int i;

	printf("Stepsize? ");
	scanf("%lf", &t);
	printf("Size of Univers? ");
	scanf("%d", &size);
	univers = malloc(size * sizeof(body));

	for (i = 0; i < size; i++) {
		printf("\n%i. Körper:\n", (i + 1));
		printf("v \n");
		scanf("%lf\n %lf\n %lf", &univers[i].v.x, &univers[i].v.y,
				&univers[i].v.z);
		printf("r \n");
		scanf("%lf\n %lf\n %lf", &univers[i].r.now.x, &univers[i].r.now.y,
				&univers[i].r.now.z);
		printf("mass ");
		scanf("%lf", &univers[i].mass);
	}
	printf("\nLength of Simulation? ");
	scanf("%lf", &endtime);
	printf("\nOutput Filename? ");
	scanf("%s", output);
//file = fopen(output, "w+");
}

FILE *init_output(char *output) {
	FILE *file;

	if ((file = fopen(output, "w+")) == NULL) {
		fprintf(stderr, "Error opening outputfile %s\n", output);
		fprintf(stderr, "%s\n", strerror(errno));
		exit(1);
	}
	fprintf(file, "{");
	return (file);
}

void end_output(FILE *file) {
	fprintf(file, "}\n");

	if (fclose(file) != 0) {
		fprintf(stderr, "Error closing output%s: \n", output);
		fprintf(stderr, "%s\n", strerror(errno));
	}
}

void printpos(double x, double y, double z) {
	fprintf(file, "{%lf,%lf,%lf}", x, y, z);
}

void printstep() {
	int i;
	fprintf(file, "{");
	printpos(univers[0].r.now.x, univers[0].r.now.y, univers[0].r.now.z);
	for (i = 1; i < size; i++) {
		fprintf(file, ",");
		printpos(univers[i].r.now.x, univers[i].r.now.y, univers[i].r.now.z);
	}
	fprintf(file, "}\n");
}

int main() {
	double time;

	user_prompt();
	file = init_output(output);

	// printout starting values
	printstep();

	setup_univers();
	printstep();

	for (time = 0; time <= endtime; time += t) {
		calc_next();
		printstep();
	}

	end_output(file);
	return 0;
}
