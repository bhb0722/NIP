#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_BUF 512
#define SPOTS 25
#define DATAS 17520

bool dispose_mem();

int interpolation(double *data);
double get_intpol(int d1, int d2, double x1, double x2);

double atod(char *str);

bool make_clp_db();
void rtrim(char *str);

/*
 * 1 : SO2
 * 2 : CO2
 * 3 : O3
 * 4 : NO2
 * 5 : PM10
*/

char spot[SPOTS][20] = { "강남구", "강동구", "강북구", "강서구", "관악구", "광진구", "구로구", "금천구", "노원구", "도봉구", "동대문구",
                      "동작구", "마포구", "서대문구", "서초구", "성동구", "성북구", "송파구", "양천구", "영등포구", "용산구", "은평구",
	                  "종로구", "중구", "중랑구"};

double **data;

int main() {
	FILE *fi;
	char fname[50];
	char check;
	int i = 0, j = 0;
	int section = 120;
	int time_delay = 1;
	double section_length;
	float t_rate = 0.9f;
	int train_index, test_index;
	double max[SPOTS] = {0., 0., 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.};
	double *avg;
	double rsq, average;
	int *count;
	
	char buf[MAX_BUF], *token;

	data = (double **)malloc(sizeof(double*) * SPOTS);
	for (int i = 0; i < SPOTS; i++) {
		data[i] = (double *)malloc(sizeof(double) * 20000);
	}

	printf("Make DB(Y?) : ");
	scanf("%c", &check);

	if (check == 'Y') {
		if (!make_clp_db()) {
			printf("error\n");
			dispose_mem();
			return 0;
		}
	}

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
				max[i] = data[i][j] > max[i] ? data[i][j] : max[i];
				token = strtok(NULL, ",");
				i++;
			}
			j++;
			i = 0;
		}

		for (i = 0; i < SPOTS; i++)
			interpolation(data[i]);

		train_index = (int)(DATAS * t_rate);
		test_index = train_index + 1;

#if 1

		for (i = 0; i < SPOTS; i++) {
			sprintf(fname, "result/%s_real_est_Ele%d_TimeDelay%d.csv", spot[i], ele, time_delay);
			FILE *fo = fopen(fname, "w");
			if (fo == NULL) {
				printf("FILE OPEN ERROR\n");
				dispose_mem();
				return 0;
			}

			avg = (double *)calloc(section, sizeof(double));
			count = (int *)calloc(section, sizeof(int));

			// Get section average
			section_length = max[i] / section;
			for (j = 0; j < train_index; j++) {
				int index = (int)(data[i][j] / section_length);
				if (index >= section)
					index = section - 1;
				avg[index] += data[i][j + time_delay];
				count[index]++;
			}

			average = 0.; rsq = 0.;

			for (j = 0; j < section; j++) {
				if (count[j] != 0) {
					avg[j] /= count[j];
				}/*
				else if(j != 0)
					avg[j] = avg[j - 1];
					*/
				fprintf(fo, "%lf\n", avg[j]);
			}

			double rmse = 0.;
			double est;

			for (j = test_index; j < DATAS; j++) {
				average += data[i][j];
			}
			average /= (DATAS - test_index);

			for (j = test_index; j < DATAS; j++) {
				int index = (int)(data[i][j - time_delay] / section_length);
				if (index >= section)
					index = section - 1;
				if (avg[index] == 0) {
					est = data[i][j - time_delay];
				}
				else {
					est = avg[index];
				}
				rmse += pow(data[i][j] - est, 2);
				rsq += pow(data[i][j] - average, 2);
				fprintf(fo, "%f,%f\n", data[i][j], est);
			}

			rsq = 1. - rmse / rsq;

			rmse /= (DATAS - test_index);
			rmse = sqrt(rmse);

			fprintf(fo, "RMSE,%lf", rmse);
			fprintf(fo, "RSQ,%lf", rsq);

			fclose(fo);

			printf("RMSE of %s : %lf / RSQ: %lf\n", spot[i], rmse, rsq);

			free(count);
			free(avg);

			sprintf(fname, "%s_cooc_ele%d_timedelay%d.csv", spot[i], ele, time_delay);
			fo = fopen(fname, "w");

			for (int fi = 0; fi < DATAS - time_delay; fi++) {
				fprintf(fo, "%f,%f\n", data[0][fi], data[0][fi + time_delay]);
			}

			fclose(fo);

		}
#endif
		fclose(fi);

	}
	
	dispose_mem();
	return 1;
}

bool dispose_mem() {
	for (int i = 0; i < SPOTS; i++) {
		free(data[i]);
	}
	free(data);
	return true;
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
		} else {
			if (f == 0) {
				ret_val *= 10.;
				ret_val += str[i] - 48;
			} else {
				ret_val += (str[i] - 48) * pow(0.1, (i - f));
			}
		}
	}

	return pos * ret_val;
}

int interpolation(double *data) {

	int i, j, k;
	int x2;
	int data_cnt;
	int cnt = 0;

	for (j = 0; j < 17520 ; j++) {
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

double get_intpol(int d1, int d2, double x1, double x2) {

	return ((d2*x1) / (d1 + d2)) + ((d1*x2) / (d1 + d2));

}

bool make_clp_db() {
	FILE *fo1, *fo2, *fo3, *fo4, *fo5, *fi;
	char fname[50];
	char buf[MAX_BUF];
	char *token;
	char date[12];
	char cur_spot[20];
	int spot_index;
	int tok_cnt;
	
	sprintf(fname, "data/1.csv");
	fo1 = fopen(fname, "w");
	sprintf(fname, "data/2.csv");
	fo2 = fopen(fname, "w");
	sprintf(fname, "data/3.csv");
	fo3 = fopen(fname, "w");
	sprintf(fname, "data/4.csv");
	fo4 = fopen(fname, "w");
	sprintf(fname, "data/5.csv");
	fo5 = fopen(fname, "w");

	strcpy(date, "2014010101\0");
	spot_index = 0;
	for (int year = 2014; year <= 2015; year++) {
		for (int i = 1; i <= 4; i++) {
			sprintf(fname, "data/%d/%d년%d분기.csv", year, year, i);
			fi = fopen(fname, "r");
			if (fi == NULL) {
				printf("[ERROR] Input file open error\n");
				fclose(fo1);
				fclose(fo2);
				fclose(fo3);
				fclose(fo4);
				fclose(fo5);
				return false;
			}
			fgets(buf, MAX_BUF, fi);
			while (fgets(buf, MAX_BUF, fi) != NULL) {
				rtrim(buf);
				token = strtok(buf, ",");

				token = strtok(NULL, ",");
				strcpy(cur_spot, token);

				token = strtok(NULL, ",");
				tok_cnt = 1;

				if (strcmp(date, token) != 0) {

					while (spot_index != SPOTS) {
						fprintf(fo1, "-999,");
						fprintf(fo2, "-999,");
						fprintf(fo3, "-999,");
						fprintf(fo4, "-999,");
						fprintf(fo5, "-999,");
						spot_index++;
					}

					strcpy(date, token);
					fprintf(fo1, "\n", token);
					fprintf(fo2, "\n", token);
					fprintf(fo3, "\n", token);
					fprintf(fo4, "\n", token);
					fprintf(fo5, "\n", token);
					spot_index = 0;
				}

				while (strcmp(cur_spot, spot[spot_index]) != 0) {
					fprintf(fo1, "-999,");
					fprintf(fo2, "-999,");
					fprintf(fo3, "-999,");
					fprintf(fo4, "-999,");
					fprintf(fo5, "-999,");
					spot_index++;
				}

				token = strtok(NULL, ",");
				while (token != NULL) {
					switch (tok_cnt) {
					case 1:
						fprintf(fo1, "%s,", token);
						break;
					case 2:
						fprintf(fo2, "%s,", token);
						break;
					case 3:
						fprintf(fo3, "%s,", token);
						break;
					case 4:
						fprintf(fo4, "%s,", token);
						break;
					case 5:
						fprintf(fo5, "%s,", token);
						break;
					}
					token = strtok(NULL, ",");
					tok_cnt++;
				}
				spot_index++;
			}
			fclose(fi);
		}
	}

	fclose(fo1);
	fclose(fo2);
	fclose(fo3);
	fclose(fo4);
	fclose(fo5);

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