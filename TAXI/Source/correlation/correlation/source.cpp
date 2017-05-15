#define _CRT_SECURE_NO_DEPRECATE
#include<iostream>
#include<cstdlib>
#include<cmath>


using namespace std;
double avg(int *point, int numberUsed);
double stdDev(int point[], int numberUsed);
double Cov(int point[], int point2[], int numberUsed, int row, int col);
void Compute(int point[], int numberUsed, int row, int col);


double cov[8][14][13] = { 0 };
double av[14][13] = { 0 };
double stdv[14][13] = { 0 };
int coun = 0;

int main() {

	FILE *fi;
	char fname[20];
	int num[3000] = {0};


	

	for (int i = 0; i <= 13; i++) {

		for (int j = 0; j <= 12; j++) {
			if (j == 12 && (i == 8 || i == 9 || i == 10)) {}
			else {
				sprintf(fname, "result/z_%d_%d.txt", j, i);
				fi = fopen(fname, "r");

				
				for (int k = 0; fscanf(fi, "%d", &num[k]) == 1; k++) {
					coun++;
				}
				av[i][j] = avg(num, coun);
				stdv[i][j] = stdDev(num, coun);

				if (i > 0 && i < 13 && j>0 && j < 12) {
						Compute(num, coun,j, i);
				}
				

				fclose(fi);
				coun = 0;
			}
			
		}
	}


	printf("average\n");
	for (int i = 13; i >= 0; i--) {
		for (int j = 0; j <= 12; j++) {
			printf("%.2f ", av[i][j]);
		}
		printf("\n");

	}

	printf("\n\n\n");
	printf("stdDev\n");

	for (int i = 13; i >= 0; i--) {
		for (int j = 0; j <= 12; j++) {
			printf("%.2f ", stdv[i][j]);
		}
		printf("\n");

	}

	printf("\n\n\n");
	printf("correlation\n");
	for (int i = 1; i <= 12; i++) {
		for (int j = 1; j <= 11; j++) {
			printf("%d, %d : ", j, i);
			for(int k =0; k <=7 ; k++)
			printf("%lf ", cov[k][i][j]);
			printf("\n");
		}
		printf("\n\n");
		

	}

	return 0;

}

double avg(int *point, int numberUsed)
{
	double sum = 0, avgValue = 0;

	for (int i = 0; i < numberUsed; i++)
	{
		sum += (double)point[i];
	}
	avgValue = sum / (double)numberUsed;
	return avgValue;  
}
	
double stdDev(int point[], int numberUsed )
{
	double avgValue;
	double sum = 0, stdValue = 0;
	avgValue = avg(point, numberUsed);
	for (int i = 0; i < numberUsed; i++)
	{
		sum += ((double)point[i] - avgValue) * ((double)point[i] - avgValue);
	}
	sum /= (double)numberUsed;
	stdValue = sqrt(sum);
	return stdValue;  
}


// correlation
double Cov(int point[], int point2[],int numberUsed,int row, int col )
{
	
	double avgValue; double avgValue2;
	double sum = 0;

	avgValue = avg(point, numberUsed);
	avgValue2 = avg(point2, numberUsed);

	for (int i = 0; i < numberUsed; i++)
	{
		sum += ((double)point[i] - avgValue) * ((double)point2[i] - avgValue2);
	}
	sum /= (double)numberUsed;

	sum /= stdDev(point, numberUsed);
	sum /= stdDev(point2, numberUsed);
	
	return sum;
}

void Compute(int point[], int numberUsed, int row, int col) {
	FILE *fp;
	char fname[20];
	int tmp[3000] = {0};
	int arr_r[8] = {-1,0,1,-1,1,-1,0,1};
	int arr_c[8] = {-1,-1,-1,0,0,1,1,1};

	for (int a = 0; a < 8; a++) {

		if (row == 11 && (col == 11 || col == 7||col == 8 || col == 9 || col == 10)) {}
		else{
			
			sprintf(fname, "result/z_%d_%d.txt", row + arr_r[a], col + arr_c[a]);
			fp = fopen(fname, "r");
			for (int k = 0; fscanf(fp, "%d", &tmp[k]) == 1; k++);
			fclose(fp);
			cov[a][col][row] = Cov(point, tmp, coun, row + arr_r[a], col + arr_c[a]);
			//printf("%lf\n", Cov(point, tmp, coun));
		}
		
	}

}