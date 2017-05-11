#define _CRT_SECURE_NO_DEPRECATE
//#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstdio>
#define SIZE 2879
#define Trate 0.7
int TR = 2016;
int TE = 863;
int dev = 25;
double v_d, v_c;
double a_d = 0, a_c = 0;

double w_d, w_c;
using namespace std;

void Cooccurrence(int *train, int *test, double *pre_val, int *ma);
void Weather(int *arr);
void DecisionTreeRegression(int* train, int *test, int *arr, double *pre_val,FILE *fi);
void Combination(int* train, double *pre_dec, double *pre_coc, FILE *fo);
void Evaluation(int *test, double *pre_dec, double *pre_coc, FILE *fo, int count);
bool readData(int col, int row, int *train, int *test, int *ma);

double tr_avg = 0.f;

int main() {
	int i, j, k,w, max;
	int *arr, *train, *test;
	FILE *fo, *fi;
	double *pre_dec, *pre_coc, *pre_out;
	double val[14][13];
	arr = (int *)calloc(SIZE, sizeof(int));
	pre_coc = (double *)calloc(SIZE, sizeof(double));
	pre_dec = (double*)calloc(SIZE, sizeof(double));
	pre_out = (double*)calloc(SIZE, sizeof(double));
	char fname[100];

	
	fi = fopen("avg.csv", "r");
	for (i = 0; i <= 13; i++)
		for (j = 0; j <= 12; j++) {
			fscanf(fi, "%lf,", &val[13 - i][j]);
		}
	fclose(fi);
	//fi = fopen("av.csv", "w");
	//TR = SIZE *Trate;		// Training rate
	//TE = SIZE - TR;
	if (Trate == 1) {
		TE = TR;
	}
	train = (int*)calloc(TR, sizeof(double));
	test = (int*)calloc(TE, sizeof(double));
	Weather(arr);
	
	for (k = 0; k <=2; k++) {		// evalution method ( 0 : combine , 1: decison tree , 2: cooccurence )
		
		if (k == 0)
			sprintf(fname, "eval_combine.csv");
		else if (k == 1)
			sprintf(fname, "eval_decision_tree.csv");
		else
			sprintf(fname, "eval_cooccurence.csv");

		fo = fopen(fname, "w");
		for (i = 13; i >=0; i--) {
			for (j = 0; j <= 12; j++) {
				


				if (!readData(i, j, train, test, &max)) {
					fprintf(fo, "0,");
					continue;
				}
				if (tr_avg < 30) {
					dev = 1;
				}
				else if (max < 100) {
					dev = 5;
				}
				else {
					dev = 25;
				}
				/*else {
					for (w = 2; w <= 53; w++) {
						if (max <= w * 50-20&& max >(w-1)*50-20) {
							dev = w;
						}
					}
				}*/
				
				DecisionTreeRegression(train, test, arr, pre_dec,fi);
				Cooccurrence(train, test, pre_coc, &max);
				Combination(train,  pre_dec, pre_coc, fo);
				Evaluation(test, pre_dec, pre_coc, fo, k);
				

				
			}
			fprintf(fo, "\n");
		}

		fclose(fo);
	}

	//fclose(fi);
	free(test);
	free(train);
	free(arr);
	free(pre_coc);
	free(pre_out);
	free(pre_dec);
	return 0;
}

bool readData(int col, int row, int *train, int *test, int *ma) {
	int val;
	char fname[20];
	int i, j = 0;
	int max = 0;
	vector<double> arr;

	sprintf(fname, "in/z_%d_%d.txt", col, row);
	FILE *fp = fopen(fname, "r");
	if (fp == NULL) {
		return false;

	}


	for (i = 0; fscanf(fp, "%d\n", &val) == 1; i++) {

		if (max < val) {
			max = val;
		}
		if (i < TR) {
			train[i] = val;
			tr_avg += val;
		}
		else {
			test[j] = val;
			j++;
		}
	}

	tr_avg /= TR;

	if (i < 2800) {
		return false;
	}

	*ma = max;

	fclose(fp);
	return true;
}

void Cooccurrence(int *train, int *test, double *pre_val, int *ma) {

	vector<int> arr;
	double *arr2;
	double *count;
	int index;
	int val;
	int len;
	int max;
	double avg = 0.0;
	char fname[20];
	int i;

	max = *ma;
	len = max / dev + 1;

	arr2 = (double*)calloc(len, sizeof(double));
	count = (double *)calloc(len, sizeof(double));

	for (i = 0; i < TR - 1; i++) {
		arr2[train[i] / dev] += train[i + 1];
		count[train[i] / dev] += 1;

	}

	for (i = 0; i <len; i++) {
		if (count[i])
			arr2[i] /= count[i];
		else if (i != 0)
			arr2[i] = arr2[i - 1]+1;
		else {
			arr2[0] = dev / 2 ;
		}
			//
	}

	for (i = 0; i < TE; i++) {
		pre_val[i] = arr2[test[i] / dev];
	}

	/* 실수 하셧네요
	
	for (i = 0; i < TR; i++) {
		pre_val[i] = arr2[train[i] / dev];
	}
	*/

	free(count);
	free(arr2);

}

void DecisionTreeRegression(int* train, int* test, int *arr, double *pre_val,FILE *fi) {

	int i, month, day, hour;
	int j = 0;

	double av = 0;
	double avg[7][24][2] = { 0 }; //day, time, rain(0,1)
	int ind[7][24][2] = { 0 };


	for (i = 0; i < TR; i++) {

		avg[i % 7][i % 24][arr[i]] += train[i];
		ind[i % 7][i % 24][arr[i]]++;
		//fprintf(fi, "%d,%d\n", train[i], arr[i]);
	}




	for (day = 0; day < 7; day++) {
		for (hour = 0; hour < 24; hour++) {
			if (ind[day][hour][1] * ind[day][hour][0] == 0) {
				double s = avg[day][hour][0] + avg[day][hour][1], t = ind[day][hour][0] + ind[day][hour][1];
				avg[day][hour][0] = s / t;
				avg[day][hour][1] = s / t;
			}
			else {

				for (i = 0; i <= 1; i++)
					avg[day][hour][i] /= ind[day][hour][i];
			}
			
		}

	}

	for (i = 0; i < TE; i++) {	//Decision Tree regression model
		j = TR + i;
		pre_val[i] = avg[i % 7][i % 24][arr[j]];
	}

}

void Weather(int *arr) {
	FILE *fp = fopen("weather.txt", "r");
	int i, month, day, hour;
	int **tmp;
	tmp = (int **)malloc(sizeof(int*) * 220);
	for (i = 0; i < 220; i++) tmp[i] = (int*)calloc(4, sizeof(int));


	for (i = 0; i < 220; i++) {
		fscanf(fp, "%d %d %d %d\n", &tmp[i][0], &tmp[i][1], &tmp[i][2], &tmp[i][3]);

	}
	int index = 0;
	int mon[5] = { 0,31,28,31,30 };

	for (month = 1; month <= 4; month++) {
		for (day = 1; day <= mon[month]; day++) {
			for (hour = 0; hour <= 23; hour++) {
				for (i = 0; i < 220; i++) {
					if (month == tmp[i][0] && day == tmp[i][1] && hour == tmp[i][2] && tmp[i][3] != 0) {
						arr[index] = 1;

					}
				}
				index++;
			}

		}
	}
	fclose(fp);
	for (i = 0; i < 220; i++) free(tmp[i]);
	free(tmp);
}

void Combination(int* train, double *pre_dec, double *pre_coc,  FILE *fo) {
	int i;
	double nq;
	double pq;
	double f_d, f_c;

	a_d = 0.f;
	a_c = 0.f;
	v_d = 0.f;
	v_c = 0.f;

	for (i = 0; i < TR; i++) {
		a_d += pre_dec[i] - train[i];
		a_c += pre_coc[i] - train[i];
	}
	a_d /= TR;
	a_c /= TR;

	for (i = 0; i < TR; i++) {
		v_d += (pre_dec[i] - train[i] - a_d)*(pre_dec[i] - train[i] - a_d);
		v_c += (pre_coc[i] - train[i] - a_c)*(pre_coc[i] - train[i] - a_c);

	}
	v_d /= TR;
	v_c /= TR;

	w_d = v_c / (v_c + v_d);
	w_c = 1.f - w_d;
}

void Evaluation(int*test, double *pre_dec, double *pre_coc, FILE *fo, int count) {

	int i;
	double val;
	double rsq, rsq2, rsq3;
	double pq = 0, nq = 0;
	double avg = 0;
	for (i = 0; i<TE; i++) {

		avg += test[i];
	}

	avg /= TE;
	for (i = 0; i < TE; i++) {
		val = w_d * pre_dec[i] + w_c * pre_coc[i];


		pq += (test[i] - val)*(test[i] - val);
		nq += (test[i] - avg)*(test[i] - avg);
	}


	rsq = 1 - pq / nq;			// combine
	pq = 0;


	for (i = 0; i < TE; i++) 
		pq += (test[i] - pre_dec[i])*(test[i] - pre_dec[i]);
	
	rsq2 = 1 - pq / nq;			// decision tree
	pq = 0;

	for (i = 0; i < TE; i++) 
		pq += (test[i] - pre_coc[i])*(test[i] - pre_coc[i]);

	rsq3 = 1 - pq / nq;			// cooccurence

								//cout << rsq << " " << rsq2 << " " << rsq3 << endl;

	if (count == 0)
		fprintf(fo, "%.4lf,", rsq);
	else if (count == 1)
		fprintf(fo, "%.4lf,", rsq2);
	else
		fprintf(fo, "%.4lf,", rsq3);
	//cout << rsq << endl;
}