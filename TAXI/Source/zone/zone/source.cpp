#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <cstdlib>
using namespace std;

int main() {
	int inum = 0;
	char fname[20];
	FILE *fp = fopen("change.txt", "r");
	int i, idx, val, pidx = -1;
	FILE *fo;
	int mat[14][13] = {0};
	int row;
	int col;


	for (i = 0; fscanf(fp, "%d", &idx) == 1; ++i) {
		row = idx % 200 - 98;
		col = idx / 200 - 70;
		
		fscanf(fp, "%d", &val);
		mat[col][row] =	1;

		if (pidx == -1 || pidx != idx) {
			sprintf(fname, "result/z_%d_%d.txt", row,col);

			if (pidx != -1)
				fclose(fo);
			fo = fopen(fname, "w");
		}

		fprintf(fo, "%d\n", val);

		pidx = idx;
	}
	printf("done\n");

	fclose(fp);
	fclose(fo);

	for (int i = 13; i >= 0; i--) {
		for (int j = 0; j <= 12; j++) {
			printf("%d ", mat[i][j]);
		}
		printf("\n");

	}


	return 0;
}