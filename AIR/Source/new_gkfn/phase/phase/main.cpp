#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "gkfn.h"

#define SELET 1
#define SPOTS 25
#define MAX_BUF 1024
#define DATAS 17520
#define DAYS 730

#define VER1 0
#define VER2 1
#if VER2
#define DAY 1
#define HOUR 0
#endif

#if 0
#define TRATE_MIN 0.7f
#define TRATE_MAX 0.9f
#define TAU_MIN 3
#define TAU_MAX 7
#define E_MIN 9
#define E_MAX 11
#define T_RATE_MIN 0.4f
#define T_RATE_MAX 0.7f
#define NOK_MIN 6
#define NOK_MAX 10
#endif

double *ts;
void calSmoothnessMeasure(int num, double threshold, int n, double* ts, int M, int &rE, int &rTau);
GKFN* predictSeries(int n, double *ts,
	int E, int Tau, int PredictionStep, double TrainingRate,
	int NumberOfKernels, int Epochs, int aq, char *spot);
#if VER1
int extract_data(char *spot_name, int ele);
int interpolation(char *file);
#endif
#if VER2
int interpolation(double *data);
bool daily_avg(double *data, double *data_days);
#endif
bool dispose_mem();
double atod(char *str);
void rtrim(char *str);
double get_intpol(int d1, int d2, double x1, double x2);

int TAU_MIN = 1;
int TAU_MAX = 10;
int E_MIN = 1;
int E_MAX = 12;
double T_RATE_MIN = 0.4;
double T_RATE_MAX = 0.6;
int NOK_MIN = 3;
int NOK_MAX = 10;
#if VER2
double **data;
double **data_days;
#endif

int main() {
	int i, j, n;
	int E, Tau, nok;
	char fname[50];
#if VER1
	char data[20];
#endif
	int year;
	GKFN *model;
	FILE *fi, *fo, *foo;
#if VER1
	ts = (double*)malloc(sizeof(double) * 20000);
#endif
	double avg;
	char spot[SPOTS][20] = { "강남구", "강동구", "강북구", "강서구", "관악구", "광진구", "구로구", "금천구", "노원구", "도봉구", "동대문구",
		"동작구", "마포구", "서대문구", "서초구", "성동구", "성북구", "송파구", "양천구", "영등포구", "용산구", "은평구",
		"종로구", "중구", "중랑구" };
	double rsq_max[39][4];

	int mode, all;
	int ele = 1;
	int ele_max;
	/*
	1 : SO2
	2 : CO
	3 : O3
	4 : NO2
	5 : PM10
	6 : PM2.5
	*/
#if VER1
	printf("!1 : All spots, 1 : select spot\n");
	scanf("%d", &all);

	if (all == 1) {
		printf("Select spot\n");
		scanf("%s", data);
	}

	printf("Select element\n");
	printf("1 : SO2\n"
		"2 : CO\n"
		"3 : O3\n"
		"4 : NO2\n"
		"5 : PM10\n"
		"6 : All of these\n");
	scanf("%d", &ele);

	printf("!1 : Preset, 1 : Range setting \n");
	scanf("%d", &mode);

	if (mode == 1) {
		printf("Training rate min(0.1~0.9) : ");
		scanf("%lf", &T_RATE_MIN);
		printf("Training rate max(0.1~0.9) : ");
		scanf("%lf", &T_RATE_MAX);

		if (T_RATE_MAX < T_RATE_MIN) {
			printf("MAX value should be bigger then MIN value\n");
			free(ts);
			return 0;
		}

		printf("Tau min : ");
		scanf("%d", &TAU_MIN);
		printf("Tau max : ");
		scanf("%d", &TAU_MAX);

		if (TAU_MAX < TAU_MIN) {
			printf("MAX value should be bigger then MIN value\n");
			free(ts);
			return 0;
		}

		printf("Embedding dimension min : ");
		scanf("%d", &E_MIN);
		printf("Embedding dimension max : ");
		scanf("%d", &E_MAX);

		if (E_MAX < E_MIN) {
			printf("MAX value should be bigger then MIN value\n");
			free(ts);
			return 0;
		}

		printf("Number of Kernals min : ");
		scanf("%d", &NOK_MIN);
		printf("Number of Kernals max : ");
		scanf("%d", &NOK_MAX);

		if (NOK_MAX < NOK_MIN) {
			printf("MAX value should be bigger then MIN value\n");
			free(ts);
			return 0;
		}

	}
	if (ele != 6) {
		i = ele;
		ele_max = ele;
	} else {
		i = 1;
		ele = 1;
		ele_max = 5;
	}

	for (i; i <= ele_max; i++) {
		for (int spot_index = 0; spot_index < 39; spot_index++) {
			if (all == 1) {
				spot_index = 100;
			} else {
				strcpy(data, spot[spot_index]);
				rsq_max[spot_index][0] = 0;
			}

			extract_data(data, i);

			sprintf(fname, "data/%s_%d.csv", data, i);
			fi = fopen(fname, "r");

			avg = 0.f;
			for (j = 1; fscanf(fi, "%lf", &ts[j]) == 1; ++j) avg += ts[j];
			fclose(fi);
			n = j - 1;
			avg /= n;

			calSmoothnessMeasure(i, 0.5f, n, ts, 1, E, Tau);

			sprintf(fname, "result/%s_%d.csv", data, i);

			fo = fopen(fname, "w");

			fprintf(fo, "E,Tau,Trate,nok,TrainRsq,TestRsq,TrainRmse,TestRmse\n");
			printf("E\tTau\tTrate\tnok\tTrainRsq\tTestRsq\t\tTrainRmse\tTestRmse\n");

			for (double trate = T_RATE_MIN; trate <= T_RATE_MAX; trate += 0.1) {
				for (nok = NOK_MIN; nok <= NOK_MAX; nok++) {
					double trrsq, trrmse, tersq, termse, ytrue, yest;
					model = predictSeries(n, ts, E, Tau, 1, trate, nok, 5, i);
					trrsq = model->getTrainRsquared();
					tersq = model->getTestRsquared();
					trrmse = model->getTrainRMSE();
					termse = model->getTestRMSE();

					// 이거
					if (isnan(trrsq) || isnan(tersq) || isnan(trrmse) || isnan(termse)) {
						delete model;
						break;
					}

					fprintf(fo, "%d,%d,%.1lf,%d,%lf,%lf,%lf,%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);
					printf("%d\t%d\t%.1lf\t%d\t%lf\t%lf\t%lf\t%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);

					if (all != 1) {
						if (rsq_max[spot_index][0] < tersq) {
							rsq_max[spot_index][0] = tersq;
							rsq_max[spot_index][1] = E;
							rsq_max[spot_index][2] = Tau;
							rsq_max[spot_index][3] = nok;
						}
					}

					delete model;
				}
			}
			fclose(fo);
		}

		if (all != 1) {
			sprintf(fname, "result/서울전체_%d.csv", i);
			foo = fopen(fname, "w");

			fprintf(foo, "측정소,E,Tau,NoK,TrainRsq\n");

			for (int ii = 0; i < 39; ii++) {
				fprintf(foo, "%s,", spot[ii]);
				fprintf(foo, "%f,", rsq_max[ii][1]);
				fprintf(foo, "%f,", rsq_max[ii][2]);
				fprintf(foo, "%f,", rsq_max[ii][3]);
				fprintf(foo, "%lf\n", rsq_max[ii][0]);
			}

			fclose(foo);
		}
	}
	free(ts);
#endif
	
#if VER2
	char buf[500];
	char *token;

	data = (double **)malloc(sizeof(double*) * SPOTS);
	for (int i = 0; i < SPOTS; i++) {
		data[i] = (double *)malloc(sizeof(double) * DATAS);
	}
#if DAY
	data_days = (double **)malloc(sizeof(double*) * SPOTS);
	for (int i = 0; i < SPOTS; i++) {
		data_days[i] = (double *)malloc(sizeof(double) * DAYS);
	}
#endif

	for (int ele = 1; ele <= 5; ele++) {
		i = 0; j = 0;
		sprintf(fname, "data/%d.csv", ele);
		fi = fopen(fname, "r");

		if (fi == NULL) {
			printf("[ERROR] file open error");
			dispose_mem();
			return 0;
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

		for (i = 0; i < SPOTS; i++)
			interpolation(data[i]);
#if DAY
		for (i = 0; i < SPOTS; i++)
			daily_avg(data[i], data_days[i]);
#endif

		for (i = 0; i < SPOTS; i++) {
#if DAY
			calSmoothnessMeasure(ele, 0.5f, DAYS, data_days[i], 1, E, Tau);

			sprintf(fname, "result/%s_%d.csv", spot[i], ele);
			fo = fopen(fname, "w");

			fprintf(fo, "E,Tau,Trate,nok,TrainRsq,TestRsq,TrainRmse,TestRmse\n");
			printf("E\tTau\tTrate\tnok\tTrainRsq\tTestRsq\t\tTrainRmse\tTestRmse\n");

			for (double trate = T_RATE_MIN; trate <= T_RATE_MAX; trate += 0.1) {
				for (nok = NOK_MIN; nok <= NOK_MAX; nok++) {
					double trrsq, trrmse, tersq, termse, ytrue, yest;
					model = predictSeries(DAYS, data_days[i], E, Tau, 1, trate, nok, 5, ele, spot[i]);
					trrsq = model->getTrainRsquared();
					tersq = model->getTestRsquared();
					trrmse = model->getTrainRMSE();
					termse = model->getTestRMSE();

					// 이거
					if (isnan(trrsq) || isnan(tersq) || isnan(trrmse) || isnan(termse)) {
						delete model;
						break;
					}

					fprintf(fo, "%d,%d,%.1lf,%d,%lf,%lf,%lf,%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);
					printf("%d\t%d\t%.1lf\t%d\t%lf\t%lf\t%lf\t%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);

					delete model;
				}
			}
#else
			calSmoothnessMeasure(ele, 0.5f, DATAS, data[i], 1, E, Tau);

			sprintf(fname, "result/%s_%d.csv", spot[i], ele);
			fo = fopen(fname, "w");

			fprintf(fo, "E,Tau,Trate,nok,TrainRsq,TestRsq,TrainRmse,TestRmse\n");
			printf("E\tTau\tTrate\tnok\tTrainRsq\tTestRsq\t\tTrainRmse\tTestRmse\n");

			for (double trate = T_RATE_MIN; trate <= T_RATE_MAX; trate += 0.1) {
				for (nok = NOK_MIN; nok <= NOK_MAX; nok++) {
					double trrsq, trrmse, tersq, termse, ytrue, yest;
					model = predictSeries(DATAS, data[i], E, Tau, 1, trate, nok, 5, i);
					trrsq = model->getTrainRsquared();
					tersq = model->getTestRsquared();
					trrmse = model->getTrainRMSE();
					termse = model->getTestRMSE();

					// 이거
					if (isnan(trrsq) || isnan(tersq) || isnan(trrmse) || isnan(termse)) {
						delete model;
						break;
					}

					fprintf(fo, "%d,%d,%.1lf,%d,%lf,%lf,%lf,%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);
					printf("%d\t%d\t%.1lf\t%d\t%lf\t%lf\t%lf\t%lf\n", E, Tau, trate, nok, trrsq, tersq, trrmse, termse);

					delete model;
				}
			}
#endif
			fclose(fo);
		}

		fclose(fi);

	}

	dispose_mem();
#endif
	return 0;
}
#if VER2

bool daily_avg(double *data, double *data_days) {
	int j = 0;
	double temp = 0.;
	for (int i = 0; i < DATAS; i++) {
		temp += data[i];
		if ((i + 1) % 24 == 23) {
			data_days[j] = (temp / 24);
			j++;
			temp = 0.;
		}
	}

	return true;
}

bool dispose_mem() {
	for (int i = 0; i < SPOTS; i++) {
		free(data[i]);
	}
	free(data);
	return true;
}

void rtrim(char *str) {
	int len = strlen(str);
	len--;
	while (str[len] == '\n' || str[len] == ' ') {
		str[len] = '\0';
		len--;
	}
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
#endif
#if VER1
int extract_data(char *spot_name, int ele) {
	FILE *fi, *fo;
	char foname[50];
	char fname[50];
	char buf[500];
	char *token;
	bool flag;

	sprintf(foname, "data/%s_%d.csv", spot_name, ele);
	fo = fopen(foname, "w");

	for (int year = 2014; year <= 2015; year++) {
		for (int i = 1; i <= 4; i++) {
			sprintf(fname, "data/%d/%d년%d분기.csv", year, year, i);
			fi = fopen(fname, "r");

			fgets(buf, 500, fi);

			while (fgets(buf, 500, fi) != NULL) {
				flag = true;
				token = strtok(buf, ",");
				for (int j = 1; j <= 9; j++) {
					token = strtok(NULL, ",");
					if (j == 1) {
						if (strcmp(token, spot_name) != 0) {
							flag = false;
							break;
						}
					}
					else {
						if (j - 2 < ele) {
							continue;
						}
						else if (j - 2 == ele) {
							if (strcmp(token, "-999") == 0)
								fprintf(fo, "0\n");
							else
								fprintf(fo, "%s\n", token);
						}
						else {
							break;
						}

					}
				}
			}

			fclose(fi);

		}
	}
	
	fclose(fo);

	interpolation(foname);
	return 1;
	

}

int interpolation(char *file) {

	FILE *fp;
	int i, j, k;
	int x2;
	int data_cnt;
	int cnt = 0;
	double data[20000];

	cnt = 0;

	printf("Adjusting '%s'\n\n", file);

	fp = fopen(file, "r");

	if (fp == NULL) {
		printf("[ERROR] Can't open file\n");
		return 0;
	}

	for (j = 0; fscanf(fp, "%lf", &data[j]) == 1; j++);

	data_cnt = j;

	printf("%d data read\n", data_cnt);

	fclose(fp);

	for (j = 0; j < data_cnt; j++) {
		if (data[j] == 0) {
			k = j;
			while (data[k] == 0) k++;
			x2 = k;
			k--;

			for (k; k >= j; k--) {
				data[k] = get_intpol(k - j + 1, x2 - k, data[j - 1], data[x2]);
				cnt++;
			}
			j = x2 - 1;

		}
	}

	printf("%d changed\n", cnt);

	fp = fopen(file, "w");

	for (j = 0; j < data_cnt; j++) {
		fprintf(fp, "%lf\n", data[j]);
	}

	fclose(fp);

	return 1;
}
#endif
#if VER2
int interpolation(double *data) {

	int i, j, k;
	int x2;
	int data_cnt;
	int cnt = 0;

	for (j = 0; j < DATAS; j++) {
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

	return 1;
}
#endif
double get_intpol(int d1, int d2, double x1, double x2) {

	return ((d2*x1) / (d1 + d2)) + ((d1*x2) / (d1 + d2));

}

GKFN* predictSeries(int n, double *ts,
	int E, int Tau, int PredictionStep, double TrainingRate,
	int NumberOfKernels, int Epochs, int aq, char *spot) {
    GKFN *model;

	model = new GKFN(n, ts, E, Tau, PredictionStep, TrainingRate, aq, spot);
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