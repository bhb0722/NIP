#include <cstdio>
#include <cstdlib>

#include "gkfn.h"

double *X, *Y;

int main() {
	int i, j, n = 7*48;
	int E, Tau, nok, ri;
	char fname[20];
	GKFN *model;
	FILE *fi, *fo;
	bool *tf = new bool[n];

	double avg;

	for (i = 0; i < n; i++)
		tf[i] = false;


	sprintf(fname, "train.txt");
	fi = fopen(fname, "r");

	Y = new double[n];
	X = new double[n];

	for (j = 1;j <=n; j++) {
		double val;

		fscanf(fi, "%lf", &val);

		while (1) {
			ri = rand() % n;

			if (!tf[ri])
				break;
		}

		tf[ri] = true;
		ri++;

		Y[ri] = val / 500.f;
		X[ri] = j;
	}
	fclose(fi);

	model = new GKFN(n, 200, 10, X, Y);


	fi = fopen("test.txt", "r");

	for (i = 1; i <= n; i++) {
		fscanf(fi, "%lf", &Y[i]);
		X[i] = i;
	}

	model->prediction_function(n, X, Y);
	
	fclose(fi);


	return 0;
}
