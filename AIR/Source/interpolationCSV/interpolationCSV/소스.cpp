#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>

#define TRUE 1
#define FALSE 0
#define NUM_FILE 6
#define MAX_DATA_SIZE 20000

double get_intpol(int d1, int d2, double x1, double x2) {

	return ((d2*x1) / (d1 + d2)) +( (d1*x2) / (d1 + d2));

}

int main(void) {

	FILE *fp;
	int i, j, k;
	int x2;
	int data_cnt;
	int cnt = 0;
	char infile[20];
	double data[MAX_DATA_SIZE];

	for (i = 1; i <= NUM_FILE; i++) {

		cnt = 0;

		if (i == 2) continue;
		
		printf("Adjusting 'aq_%d'\n\n", i);

		sprintf(infile, "aq_%d.csv", i);

		fp = fopen(infile, "r");

		if (fp == NULL) {
			printf("[ERROR] Can't open file\n");
			return FALSE;
		}

		for (j = 0; fscanf(fp, "%lf", &data[j]) == 1; j++);

		data_cnt = j;

		printf("%d data readed\n", data_cnt);

		fclose(fp);

		for (j = 0; j < data_cnt; j++) {
			if (data[j] == 0) {
				k = j;
				while (data[k] == 0) k++;
				x2 = k;
				k--;

				for (k; k >= j; k--) {
					data[k] = get_intpol(k - j+1, x2 - k, data[j-1], data[x2]);
					//printf("changed\n");
					cnt++;
				}
				j = x2 - 1;

			}
		}

		printf("%d changed\n", cnt);

		fp = fopen(infile, "w");

		for (j = 0; j < data_cnt; j++) {
			fprintf(fp, "%lf\n", data[j]);
		}

		fclose(fp);
		
	}

}