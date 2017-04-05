#include <cstdio>
#include <cstdlib>

#include "gkfn.h"

#define HOUR 1
#define DAY 0
#define SELET 1

#define TRATE_MIN 0.7f
#define TRATE_MAX 0.9f
#define TAU_MIN 3
#define TAU_MAX 9
#define E_MIN 9
#define E_MAX 11
#define T_RATE_MIN 0.7f
#define T_RATE_MAX 0.9f
#define NOK_MIN 6
#define NOK_MAX 10

double *ts;
void calSmoothnessMeasure(int num, double threshold, int n, double* ts, int M, int &rE, int &rTau);
GKFN* predictSeries(int n, double *ts,
	int E, int Tau, int PredictionStep, double TrainingRate,
	int NumberOfKernels, int Epochs, int aq);

int main() {
	int i, j, n;
	int E, Tau, nok;
	char fname[20];
	GKFN *model;
	FILE *fi, *fo;
	ts = (double*)malloc(sizeof(double) * 20000);
	double avg;


	
	/*
	1 : pm 1.0
	2 : pm 2.5
	3 : O3
	4 : NO2
	5 : NO
	6 : SO2
	*/

	for (i = 1; i <=6; i++) { 

		if (i == 2) continue;

#if HOUR
		//sprintf(fname, "data/lotto.txt");
		sprintf(fname, "data/aq_%d.csv", i);
#endif
#if DAY
		sprintf(fname, "data/aqd_%d.csv", i);
#endif
		fi = fopen(fname, "r");
		

		avg = 0.f;
		for (j = 1; fscanf(fi, "%lf", &ts[j]) == 1; ++j) avg += ts[j];
		fclose(fi);
		n = j - 1;
		avg /= n;
		
		//if (avg < 2.f) continue;
#if SELET
		calSmoothnessMeasure(i,0.5f, n, ts, 1, E, Tau);
#else
		E = 6;
		Tau = 1;
#endif
#if HOUR
		sprintf(fname, "result/z_%d.csv", i);
#endif
#if DAY
		sprintf(fname, "result/zd_%d.csv", i);
#endif
		fo = fopen(fname, "w");
		
		fprintf(fo, "E,Tau,Trate,nok,TrainRsq,TestRsq,TrainRmse,TestRmse\n");
		printf("E,Tau,Trate,nok,TrainRsq,TestRsq,TrainRmse,TestRmse\n");

		for (double trate = T_RATE_MIN; trate <= T_RATE_MAX; trate += 0.1) {
			for (nok = NOK_MIN; nok <= NOK_MAX; nok++) {
				double trrsq, trrmse, tersq, termse, ytrue, yest;
				model = predictSeries(n, ts, E, Tau, 1, trate, nok, 5, i);
				trrsq = model->getTrainRsquared();
				tersq = model->getTestRsquared();
				trrmse = model->getTrainRMSE();
				termse = model->getTestRMSE();

				// ÀÌ°Å
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
	int NumberOfKernels, int Epochs, int aq) {
    GKFN *model;

	model = new GKFN(n, ts, E, Tau, PredictionStep, TrainingRate, aq);
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

#if HOUR
	sprintf(fname, "sm/sm_%d.csv", num);
#endif
#if DAY
	sprintf(fname, "sm/smd_%d.csv", num);
#endif

	ifp = fopen(fname, "w");

	fprintf(ifp, "E,Tau,Smoothness\n");


	for (E = E_MIN; E <= E_MAX; E++) {
		for (Tau = TAU_MIN; Tau <= TAU_MAX; Tau++) {
			sm = 0.f;
			Ns = n - (E - 1) * Tau - 1;

			if (Ns < 100)
				continue;


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
			printf( "%d,%d,%lf\n", E, Tau, sm);
			

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