#include <cstdio>
#include <cstdlib>

#include "gkfn.h"

double *ts;
void calSmoothnessMeasure(int n, double* ts, int M);

int main() {
	int i, j, n;
	int E, Tau;
	char fname[20];
	FILE *fi, *fo;
	GKFN *model;
	ts = (double*)malloc(sizeof(double) * 10000);

	// 0.9 : 1¹ø
	// 0.4~0.5 : 3¹ø
	// 0. : 5¹ø
	
	for (i = 3; i <=3; i++) { 
		sprintf(fname, "data/z_%d.txt", i);
		fi = fopen(fname, "r");
		

		for (j = 1; fscanf(fi, "%lf", &ts[j]) == 1; ++j);
		n = j - 1;

		//calSmoothnessMeasure(n, ts, 1);

		
		sprintf(fname, "result/z_%d.csv", i);
		fo = fopen(fname, "w");
		
		for (E = 10; E <= 10; E++) {
			for (Tau = 1; Tau <= 4; Tau++) {
				model = new GKFN(n,ts, E, Tau, 1, 0.8);
				model->learn(10, 5);

				printf("----------\n");
				printf("E = %d, Tau = %d\n", E, Tau);
				printf("Rsquared = %lf\n", model->getRsquared());
				printf("RMSE = %lf\n", model->getRMSE());
				printf("\n");

				fprintf(fo,"%d,%d,%lf,%lf\n", E,Tau,model->getRMSE(), model->getRsquared());

				delete model;
			}
		}

		for (E = 24; E <= 24; E++) {
			for (Tau = 1; Tau <= 1; Tau++) {
				model = new GKFN(n, ts, E, Tau, 1, 0.8);
				model->learn(10, 1);

				printf("----------\n");
				printf("E = %d, Tau = %d\n", E, Tau);
				printf("Rsquared = %lf\n", model->getRsquared());
				printf("RMSE = %lf\n", model->getRMSE());
				printf("\n");

				fprintf(fo, "%d,%d,%lf,%lf\n", E, Tau, model->getRMSE(), model->getRsquared());

				delete model;
			}
		}
		
		fclose(fo);
		

		fclose(fi);
		
	}

	free(ts);
	
	return 0;
}

void calSmoothnessMeasure(int n, double* ts, int M) {
	int i, j, k;
	FILE *ifp;
	char ch;
	int E, Tau, Ns;
	double **x, *y;
	int l;
	double dist, mdist;
	double sm, msm;
	int mE = 0, mTau;

	ifp = fopen("SmoothnessMeasure.csv", "w");

	fprintf(ifp, "E,Tau,Smoothness\n");

	for (E = 2; E <= 10; E++) {
		for (Tau = 1; Tau <= 24; Tau++) {
			sm = 0.f;
			Ns = n - (E - 1) * Tau - 1;

			x = (double**)calloc(Ns, sizeof(double*));
			y = (double*)calloc(Ns, sizeof(double));


			k = 1 + (E - 1) * Tau;

			for (i = 1; i <= Ns; i++)
			{
				x[i] = (double*)calloc(E, sizeof(double));
				y[i] = ts[k + M];

				for (j = 0; j < E; j++) {
					x[i][j] = ts[k - j * Tau];
				}

				k++;
			}

			for (i = 1; i <= Ns; i++) {
				mdist = -1.f;
				for (j = 1; j <= Ns; j++) {
					if (i == j) continue;

					dist = 0.f;

					for (k = 0; k < E; k++)
						dist += (x[i][k] - x[j][k])*(x[i][k] - x[j][k]);

					if (dist != 0.f && (mdist == -1.f || mdist > dist)) {
						mdist = dist;
						l = j;
					}
				}

				sm += abs(y[i] - y[l]) / sqrt(mdist);
			}

			sm = 1.f - sm / Ns;
			fprintf(ifp, "%d,%d,%lf\n", E, Tau, sm);
			printf("%d,%d,%lf\n", E, Tau, sm);

			if (mE == 0 || msm < sm) {
				msm = sm;
				mE = E;
				mTau = Tau;
			}
		}
	}

	fclose(ifp);
}