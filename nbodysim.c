/*
 * nbodysim.c
 *
 * An fast n-body simulation 
 *
 *  Created on: 09.05.2014
 *      Author: robin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Global Vars
double t;
double gc = 6.67384 * 10e-11;
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
	matrix a;
	vector v;
	matrix r;
	double mass;
	matrix f;
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
		mgrav = gc * univers[i].mass * univers[body].mass;

		tmp = mgrav * (dx / diff);
		univers[body].f.now.x += tmp;
		univers[i].f.now.x += -1 * tmp;
		tmp = mgrav * (dy / diff);
		univers[body].f.now.y += tmp;
		univers[i].f.now.y += -1 * tmp;
		tmp = mgrav * (dz / diff);
		univers[body].f.now.z += tmp;
		univers[i].f.now.z += -1 * tmp;
	}

	printf("%lf", univers[body].f.now.x);
	printf("%lf", univers[body].f.now.y);
	printf("%lf", univers[body].f.now.z);
}

void leap(int body) {
	univers[body].r.new.x = 2 * univers[body].r.now.x - univers[body].r.old.x
			+ (univers[body].f.now.x / univers[body].mass) * t * t;
	univers[body].r.new.y = 2 * univers[body].r.now.y - univers[body].r.old.y
			+ (univers[body].f.now.y / univers[body].mass) * t * t;
	univers[body].r.new.z = 2 * univers[body].r.now.z - univers[body].r.old.z
			+ (univers[body].f.now.z / univers[body].mass) * t * t;
}

void setup_next() {
	int i;
	for (i = 0; i < size; i++) {
		univers[i].f.now.x = 0;
		univers[i].f.now.y = 0;
		univers[i].f.now.z = 0;
		univers[i].r.old.x = univers[i].r.now.x;
		univers[i].r.old.y = univers[i].r.now.y;
		univers[i].r.old.z = univers[i].r.now.z;
		univers[i].r.now.x = univers[i].r.new.x;
		univers[i].r.now.y = univers[i].r.new.y;
		univers[i].r.now.z = univers[i].r.new.z;
	}
}

void setup_univers() {
	int i;

	for (i = 0; i < size; i++) {
		univers[i].f.now.x = 0;
		univers[i].f.now.y = 0;
		univers[i].f.now.z = 0;
	}

	for (i = 0; i < size; i++) {
		sum_grav(i);
	}

	for (i = 0; i < size; i++) {
		univers[i].r.new.x = univers[i].r.now.x + univers[i].v.x * t
				+ 0.5 * univers[i].a.now.x * t * t;
		univers[i].r.new.y = univers[i].r.now.y + univers[i].v.y * t
				+ 0.5 * univers[i].a.now.y * t * t;
		univers[i].r.new.z = univers[i].r.now.z + univers[i].v.z * t
				+ 0.5 * univers[i].a.now.z * t * t;
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
		printf("\n");
		printf("mass ");
		scanf("%lf", &univers[i].mass);
	}
	printf("\nLength of Simulation? ");
	scanf("%lf", &endtime);
	printf("\nOutput Filename? ");
	scanf("%s", output);
	file = fopen(output, "w+");
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
	fprintf(file, "}");
}
int main() {
	double time;

	user_prompt();

	setup_univers();
	fprintf(file, "{");
	printstep();
	for (time = 0; time <= endtime; time += t) {
		calc_next();
		fprintf(file, ",");
		printstep();
	}
	fprintf(file, "}");
	return 0;
}
