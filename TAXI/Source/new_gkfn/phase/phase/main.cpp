#include <cstdio>
#include <cstdlib>

#include "gkfn.h"

double *ts;
void calSmoothnessMeasure(int num, double threshold, int n, double* ts, int M, int &rE, int &rTau);
GKFN* predictSeries(int n, double *ts,
	int E, int Tau, int PredictionStep, double TrainingRate,
	int NumberOfKernels, int Epochs);

int main() {
	int i, j, n;
	int E, Tau, nok;
	char fname[20];
	GKFN *model;
	FILE *fi, *fo;
	ts = (double*)malloc(sizeof(double) * 10000);
	double avg;


	// 0.9 : 1번
	// 0.4~0.5 : 3번
	// 0. : 5번
	
	for (i = 104; i <=179; i++) { 
		sprintf(fname, "data/z_%d.txt", i);
		fi = fopen(fname, "r");
		

		avg = 0.f;
		for (j = 0; fscanf(fi, "%lf", &ts[j]) == 1; +j) avg += ts[j];
		fclose(fi);
		n = j - 1;
		avg /= n;
		
		if (avg < 2.f) continue;

		calSmoothnessMeasure(i,0.4f, n, ts, 1, E, Tau);

		/*E = 10;
		Tau = 3;*/


		sprintf(fname, "result/z_%d.csv", i);
		fo = fopen(fname, "w");
		
		fprintf(fo, "E,Tau,Trate,nok,TrainRsq,TestRsq,TrainRmse,TestRmse\n");
		printf("E,Tau,Trate,nok,TrainRsq,TestRsq,TrainRmse,TestRmse\n");

		for (double trate = 0.8; trate <= 0.8; trate += 0.1) {
			for (nok = 7; nok <= 13; nok++) {
				double trrsq, trrmse, tersq, termse;
				model = predictSeries(n, ts, E, Tau, 1, trate, nok, 5);
				trrsq = model->getTrainRsquared();
				tersq = model->getTestRsquared();
				trrmse = model->getTrainRMSE();
				termse = model->getTestRMSE();

				// 이거
				if (isnan(trrsq) || isnan(tersq) || isnan(trrmse) || isnan(termse)) {
					delete model;
					break;
				}

				fprintf(fo, "%d,%d,%.1lf,%d,%lf,%lf,%lf,%lf\n",E,Tau,trate,nok,trrsq,tersq,trrmse,termse);
				printf( "%d,%d,%.1lf,%d,%lf,%lf,%lf,%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);

				delete model;
			}
		}

		fclose(fo);
		
		
	}

	free(ts);
	
	return 0;
}

GKFN* predictSeries(int n, double *ts,
	int E, int Tau, int PredictionStep, double TrainingRate,
	int NumberOfKernels, int Epochs) {
    GKFN *model;

	model = new GKFN(n, ts, E, Tau, PredictionStep, TrainingRate);
	model->learn(NumberOfKernels, Epochs);

	return model;
}

void calSmoothnessMeasure(int num, double threshold, int n, double* ts, int M, int &rE, int &rTau) {
	int i, j, k;
	FILE *ifp;
	char ch, fname[20];
	int E, Tau, Ns;
	double **x, *y;
	int l;
	double dist, mdist;
	double sm, msm;
	int mE = 0, mTau;

	rTau = 0;
	sprintf(fname, "sm/sm_%d.csv", num);
	ifp = fopen(fname, "w");

	fprintf(ifp, "E,Tau,Smoothness\n");


	for (E = 2; E <= 10; E++) {
		for (Tau = 1; Tau <= 24; Tau++) {
			sm = 0.f;
			Ns = n - (E - 1) * Tau - 1;

			x = (double**)calloc(Ns+1, sizeof(double*));
			y = (double*)calloc(Ns+1, sizeof(double));

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
			

			for (i = 1; i <= Ns; i++)
				free(x[i]);
			free(x);
			free(y);

			if (rTau == 0 && sm > threshold) {
				rE = E;
				rTau = Tau;
				printf("%d,%d,%lf\n", E, Tau, sm);
				return;
			}

			if (mE == 0 || msm < sm) {
				msm = sm;
				mE = E;
				mTau = Tau;
			}
		}
	}

	rE = mE;
	rTau = mTau;
	printf("%d,%d,%lf\n", mE, mTau, msm);

	fclose(ifp);
}