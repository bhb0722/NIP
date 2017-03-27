#include <cstdio>
#include <cstdlib>

#include "gkfn.h"

int main() {
	int i, j, M, k, n = 0;
	FILE *ifp;
	char ch;

	int E, Tau, Ns;

	double ts[10000];
	double **x, *y;

	int l;
	double dist, mdist;
	double sm, msm;
	int mE=0, mTau;


	/*printf("prediction time = ");
	scanf("%d", &M);

	ifp = fopen("sample2.txt", "r");

	for (i = 0; (ch = getc(ifp)) != EOF; ++i) {
		fscanf(ifp, "%lf", &ts[i]);
		n++;
	}
	fclose(ifp);

	ifp = fopen("output2.txt", "w");


	for (E = 2; E <= 40; E++) {
		for (Tau = 1; Tau <= 24; Tau++) {
			sm = 0.f;
			Ns = n - (E - 1) * Tau - 1;

			x = (double**)calloc(Ns, sizeof(double*));
			y = (double*)calloc(Ns, sizeof(double));


			k = (E - 1) * Tau;

			for (i = 0; i < Ns; i++)
			{
				x[i] = (double*)calloc(E, sizeof(double));
				y[i] = ts[k + M];

				for (j = 0; j < E; j++) {
					x[i][j] = ts[k - j * Tau];
				}

				k++;
			}

			for (i = 0; i < Ns; i++) {
				mdist =- 1.f;
				for (j = 0; j < Ns; j++) {
					if (i == j) continue;

					dist = 0.f;

					for (k = 0; k < E; k++)
						dist += (x[i][k] - x[j][k])*(x[i][k] - x[j][k]);

					if (mdist == -1.f || mdist > dist) {
						mdist = dist;
						l = j;
					}
				}

				if (mdist == 0) break;
				sm += abs(y[i] - y[l]) / sqrt(mdist);
			}

			if (i != Ns)
				continue;

			sm = 1.f - sm / Ns;
			fprintf(ifp, "E = %d   Tau = %d   Smoothness = %.4lf\n", E, Tau, sm);

			if (mE == 0 || msm < sm) {
				msm = sm;
				mE = E;
				mTau = Tau;

				
			}
		}
	}

	fclose(ifp);*/

	GKFN *model = new GKFN("input.txt", 6, 3, 1, 0.8);
	model->learn(20,20);
	printf("Rsquared = %lf\n", model->getRsquared());
	printf("RMSE = %lf\n", model->getRMSE());
	delete model;

	return 0;
}