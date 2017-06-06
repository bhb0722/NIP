#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <string.h>
#include "gkfn.h"

#define DEBUG 0
// If DAY is 1, learn & predict daily, else, by an hour
#define DAY 1

#define SPOTS 25
#define N_DATA 17520
#define N_DAYS 730

#define MAX_BUF 1024

int E_MIN = 1, E_MAX = 15;
int TAU_MIN = 1, TAU_MAX = 7;	
int NOK_MIN = 5, NOK_MAX = 15;
float SMTHRESHOLD = 0.5f;
double TRATE_MIN = 0.8, TRATE_MAX = 0.8;

double **data;
double *rsq_max;
#if DAY
double **daily_data;
#endif

GKFN* predictSeries(int n, double *ts,
	int E, int Tau, int PredictionStep, double TrainingRate,
	int NumberOfKernels, int Epochs, int aq, char *spot);
void calSmoothnessMeasure(int num, double threshold, int n, double* ts, int M, int &rE, int &rTau);
void readData(int element);
void interpolation(double *data);
double get_intpol(int d1, int d2, double x1, double x2);
void disposeMem();
double atod(char *str);
void rtrim(char *str);
void scaling(double *data, int num);
#if DAY
void getDailyAvg(int spot);
#endif

int main() {

	int E, Tau, nok;
	char elements[5][10] = { "SO2", "CO", "O3", "NO2", "PM10" };
	char spots[SPOTS][20] = { "강남구", "강동구", "강북구", "강서구", "관악구", "광진구", "구로구", "금천구", "노원구", "도봉구", "동대문구",
		"동작구", "마포구", "서대문구", "서초구", "성동구", "성북구", "송파구", "양천구", "영등포구", "용산구", "은평구",
		"종로구", "중구", "중랑구" };
	char fname[50];
	FILE *fo;
	GKFN *model;

	// Initialize Array
	data = (double **)calloc(SPOTS, sizeof(double *));
	for (int i = 0; i < SPOTS; i++) {
		data[i] = (double *)calloc(N_DATA+1, sizeof(double));
	}
	rsq_max = (double *)calloc(SPOTS, sizeof(double));
#if DAY
	daily_data = (double **)calloc(SPOTS, sizeof(double *));
	for (int i = 0; i < SPOTS; i++) {
		daily_data[i] = (double *)calloc(N_DAYS+1, sizeof(double));
	}
#endif
	// For each element
	for (int element = 3; element <= 5; element++) {

		readData(element);

		for (int spotIndex = 0; spotIndex < SPOTS; spotIndex++) {

			printf("%s %s\n", spots[spotIndex], elements[element-1]);
			interpolation(data[spotIndex]);

#if DAY
			getDailyAvg(spotIndex);
			scaling(daily_data[spotIndex], N_DAYS);
			calSmoothnessMeasure(element, SMTHRESHOLD, N_DAYS, daily_data[spotIndex], 1, E, Tau);
#else
			scaling(data[spotIndex], N_DATA);
			calSmoothnessMeasure(element, SMTHRESHOLD, N_DATA, data[spotIndex], 1, E, Tau);
#endif

			sprintf(fname, "result/%s_%s.csv", spots[spotIndex], elements[element-1]);
			fo = fopen(fname, "w");

			fprintf(fo, "E,Tau,Trate,nok,TrainRsq,TestRsq,TrainRmse,TestRmse\n");
			printf("E\tTau\tTrate\tnok\tTrainRsq\tTestRsq\t\tTrainRmse\tTestRmse\n");

			for (double trate = TRATE_MIN; trate <= TRATE_MAX; trate += 0.1) {
				for (nok = NOK_MIN; nok <= NOK_MAX; nok++) {
					double trrsq, trrmse, tersq, termse, ytrue, yest;
#if DAY
					model = predictSeries(N_DAYS, daily_data[spotIndex], E, Tau, 1, trate, nok, 5, element, spots[spotIndex]);
#else
					model = predictSeries(N_DATA, data[spotIndex], E, Tau, 1, trate, nok, 5, element, spots[spotIndex]);
#endif
					trrsq = model->getTrainRsquared();
					tersq = model->getTestRsquared();
					trrmse = model->getTrainRMSE();
					termse = model->getTestRMSE();

					rsq_max[spotIndex] = rsq_max[spotIndex] > tersq ? rsq_max[spotIndex] : tersq;

					if (isnan(trrsq) || isnan(tersq) || isnan(trrmse) || isnan(termse)) {
						delete model;
						break;
					}

					fprintf(fo, "%d,%d,%.1lf,%d,%lf,%lf,%lf,%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);
					printf("%d\t%d\t%.1lf\t%d\t%lf\t%lf\t%lf\t%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);

					delete model;
				}
			}

			fclose(fo);

		}
		double avg_rsq = 0.;

		for (int i = 0; i < SPOTS; i++) {
			avg_rsq += rsq_max[i];
		}
		avg_rsq /= SPOTS;

		printf("Average Rsq is : %lf\n", avg_rsq);

	}
	disposeMem();

}

void interpolation(double *data) {

	int i, j, k;
	int x2;
	int data_cnt;
	int cnt = 0;
	double max = 0., min = 9999.;
	double scale;

	for (j = 1; j <= N_DATA; j++) {
		if (data[j] == 0.) {
			k = j;
			while (data[k] == 0.) k++;
			x2 = k;
			k--;

			for (k; k >= j; k--) {
				data[k] = get_intpol(k - j + 1, x2 - k, data[j - 1], data[x2]);
				cnt++;
			}
			j = x2 - 1;

		}
	}

	//printf("%d changed\n", cnt);
}

double get_intpol(int d1, int d2, double x1, double x2) {

	return ((d2*x1) / (d1 + d2)) + ((d1*x2) / (d1 + d2));

}

void scaling(double *data, int num) {
	double max = 0., min = 9999.;
	double scale;

	for (int i = 1; i <= num; i++) {
		max = max >= data[i] ? max : data[i];
		min = min <= data[i] ? min : data[i];
	}

	scale = max - min;

	for (int i = 1; i <= num; i++) {
		data[i] -= min;
		data[i] /= scale;
	}
}

void readData(int element) {
	char fname[50], buf[MAX_BUF];
	char *token;
	int i = 0, j = 1;
	double max[SPOTS], min[SPOTS];
	double scale;
	FILE *fi;

	for (int k = 0; k < SPOTS; k++) {
		max[k] = 0.;
		min[k] = 9999.;
	}

	sprintf(fname, "data/%d.csv", element);
	fi = fopen(fname, "r");

	if (fi == NULL) {
		printf("[ERROR] file open error");
		disposeMem();
		return ;
	}

	while (fgets(buf, MAX_BUF, fi) != NULL) {
		rtrim(buf);
		token = strtok(buf, ",");
		while (token != NULL) {
			data[i][j] = atod(token);
			token = strtok(NULL, ",");
			i++;
		}
		j++;
		i = 0;
	}

}

#if DAY
void getDailyAvg(int spot) {
	double temp = 0.;
	int j = 1;
	
	for (int i = 1; i <= N_DATA; i++) {
		temp += data[spot][i];
		if (i % 24 == 0) {
			daily_data[spot][j] = temp / 24;
			j++;
			temp = 0.;
		}
	}

}
#endif

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


	for (E = E_MIN; E <= E_MAX; E++) {
		for (Tau = TAU_MIN; Tau <= TAU_MAX; Tau++) {
			sm = 0.f;
			Ns = n - (E - 1) * Tau - 1;

			if (Ns < 100)
				continue;


			x = (double**)calloc(Ns + 1, sizeof(double*));
			y = (double*)calloc(Ns + 1, sizeof(double));

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

GKFN* predictSeries(int n, double *ts, int E, int Tau, int PredictionStep, double TrainingRate,
	int NumberOfKernels, int Epochs, int aq, char *spot) {
	GKFN *model;

	model = new GKFN(n, ts, E, Tau, PredictionStep, TrainingRate, aq, spot);
	model->learn(NumberOfKernels, Epochs);

	return model;
}

void disposeMem() {
	for (int i = 0; i < SPOTS; i++)	free(data[i]);
	free(data);
#if DAY
	for (int i = 0; i < SPOTS; i++)	free(daily_data[i]);
	free(daily_data);
#endif
}

double atod(char *str) {
	int len = strlen(str);
	int i;
	double ret_val = 0., pos = 1.;
	int f = 0;

	if (str[0] == '-')
		return 0.;

	for (i = 0; i < len; i++) {
		if (str[i] == '.') {
			f = i;
		}
		else {
			if (f == 0) {
				ret_val *= 10.;
				ret_val += str[i] - 48;
			}
			else {
				ret_val += (str[i] - 48) * pow(0.1, (i - f));
			}
		}
	}

	return pos * ret_val;
}

void rtrim(char *str) {
	int len = strlen(str);
	len--;
	while (str[len] == '\n' || str[len] == ' ') {
		str[len] = '\0';
		len--;
	}
}