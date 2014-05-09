/*
 * nbodysim.c
 *
 *  Created on: 09.05.2014
 *      Author: robin
 */

#include <stdio.h>
#include <math.h>

//Global Vars
double t;
double y = 6.67384 * 10^(-11);
int size;
double time;
double endtime;
char output[100];
const FILE *file = fopen(output, "a+");

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

	for (i = body; i < size; i++) {
		double diff = sqrt((univers[i].r.now.x - univers[body].r.now.x)^2 + (univers[i].r.now.y - univers[body].r.now.y)^2 + (univers[i].r.now.z - univers[body].r.now.z) ^2)^3;

		double tmp = (univers[i].r.now.x - univers[body].r.now.x) / diff;
		univers[body].f.now.x += tmp;
		univers[i].f.now.x += -1 * tmp;
		tmp = (univers[i].r.now.y - univers[body].r.now.y) / diff;
		univers[body].f.now.y += tmp;
		univers[i].f.now.y += -1 * tmp;
		tmp = (univers[i].r.now.z - univers[body].r.now.z) / diff;
		univers[body].f.now.z += tmp;
		univers[i].f.now.z += -1 * tmp;
	}

	univers[body].f.now.x *= y * univers[body].mass;
	univers[body].f.now.y *= y * univers[body].mass;
	univers[body].f.now.z *= y * univers[body].mass;
}

void a() {
	int i;

	for (i = 0; i < size; i++) {
		univers[i].a.now.x = univers[i].f.now.x * (1 / univers[i].mass);
		univers[i].a.now.y = univers[i].f.now.y * (1 / univers[i].mass);
		univers[i].a.now.z = univers[i].f.now.z * (1 / univers[i].mass);
	}
}

void leap(int body) {
	univers[body].r.new.x = 2 * univers[body].r.now.x - univers[body].r.old.x
			+ univers[body].a.now.x * t * t;
	univers[body].r.new.y = 2 * univers[body].r.now.y - univers[body].r.old.y
			+ univers[body].a.now.y * t * t;
	univers[body].r.new.z = 2 * univers[body].r.now.z - univers[body].r.old.z
			+ univers[body].a.now.z * t * t;
}

void setup_next() {
	int i;
	for (i = 0; i < size; i++) {
		univers[i].r.new.x = univers[i].r.now.x + univers[i].v.x * t
				+ 0.5 * univers[i].a.now.x * (t*t);
		univers[i].r.new.y = univers[i].r.now.y + univers[i].v.y * t
				+ 0.5 * univers[i].a.now.y * (t*t);
		univers[i].r.new.z = univers[i].r.now.z + univers[i].v.z * t
				+ 0.5 * univers[i].a.now.z * (t*t);
		univers[i].a.now.x = univers[i].a.old.x;
		univers[i].a.now.y = univers[i].a.old.y;
		univers[i].a.now.z = univers[i].a.old.z;
		univers[i].a.new.x = univers[i].a.now.x;
		univers[i].a.new.y = univers[i].a.now.y;
		univers[i].a.new.z = univers[i].a.now.z;
	}
}

void setup_univers() {
	int i;

	for (i = 0; i < size; i++) {
		sum_grav(i);
	}

	a();

	for (i = 0; i < size; i++) {
		univers[i].r.new.x = univers[i].r.now.x + univers[i].v.x * t
				+ 0.5 * univers[i].a.now.x * (t ^ 2);
		univers[i].r.new.y = univers[i].r.now.y + univers[i].v.y * t
				+ 0.5 * univers[i].a.now.y * (t ^ 2);
		univers[i].r.new.z = univers[i].r.now.z + univers[i].v.z * t
				+ 0.5 * univers[i].a.now.z * (t ^ 2);
		univers[i].a.now.x = univers[i].a.old.x;
		univers[i].a.now.y = univers[i].a.old.y;
		univers[i].a.now.z = univers[i].a.old.z;
		univers[i].a.new.x = univers[i].a.now.x;
		univers[i].a.new.y = univers[i].a.now.y;
		univers[i].a.new.z = univers[i].a.now.z;
	}
}

void calc_next() {
	int i;

	setup_next();
int i;
	for (i = 0; i < size; i++) {
		sum_grav(body);
	}

	a();

	for (i = 0; i < size; i++) {
		leap(i);
	}
}

void user_prompt() {

	int i;

	printf("Stepsize? ");
	scanf("%lf", &t);
	printf("\n Size of Univers?");
	scanf("%d", &size);
	univers = malloc(size * sizeof(body));

	for (i = 0; i < size; i++) {
		printf("\n%i. Körper:\n", (i - 1));
		printf("v ");
		scanf("%lf %lf %lf", &univers[i].v.x, &univers[i].v.y, &univers[i].v.z);
		printf("\n");
		printf("r ");
		scanf("%lf %lf %lf", &univers[i].r.now.x, &univers[i].r.now.y,
				&univers[i].r.now.z);
		printf("\n");
		printf("mass ");
		scanf("%lf", &univers[i].mass);
	}
	printf("\nTime of Simulation?");
	scanf("%lf", &endtime);
	printf("\nOutput Filename?");
	fgets(output, 100, stdin);
}

void main() {
	double i;

	user_prompt();

	setup_univers();
	for (i = 0; i <= t; i = i + t) {
		calc_next();
		fprintf(file, "{");
		for (i = 0; i < size; i++) {
			fprintf(file, "{");
			double x = univers[i].r.now.x;
			double y = univers[i].r.now.y;
			double z = univers[i].r.now.z;
			fprintf(file, "%lf,", x);
			fprintf(file, "%lf,", y);
			fprintf(file, "%lf", z);
			fprintf(file, "},");
		}
		fprintf(file, "}");
		if( i == (size-1)) {
		fprintf(file, ",");
		}
	}
}
