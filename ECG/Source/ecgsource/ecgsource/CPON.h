#ifndef CPON_H
#define CPON_H

#include "cpon/betadist.h"
#include "cpon/GlobalFunction.h"
#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;

#define PRINT_CENT 0 // 금방됨
#define PRINT_KERNEL 0 // 금방됨
#define PRINT_HIST 0 // 조금걸림
#define PRINT_PARAM 0 // 걸림
#define PRINT_UNCERTAIN_SAMPLE 1 //

#pragma warning(disable:4996)
class CPON {


public:
	//Instructors
	CPON() {}
	CPON(int cpos, int cneg, double alpha) {
		this->cpos = cpos;
		this->cneg = cneg;
		this->alpha = alpha;
	}
	void clear() {
		int i, j;

		if (!is_train) return;

		for (i = 0; i < npos; i++)
			free(posk[i]);

		for (i = 0; i < nneg; i++)
			free(negk[i]);

		for (i = 0; i < cpos; i++)
		{
			free(posm[i]);

			for (j = 0; j < DIM; j++) {
				free(posv[i][j]);
				free(pos_incov[i][j]);
			}

			free(posv[i]);
			free(pos_incov[i]);
		}

		for (i = 0; i < cneg; i++) {
			free(negm[i]);

			for (j = 0; j < DIM; j++) {
				free(negv[i][j]);
				free(neg_incov[i][j]);
			}

			free(negv[i]);
			free(neg_incov[i]);
		}

		free(posk);
		free(negk);
		free(posm);
		free(posv);
		free(negm);
		free(negv);
		free(posp);
		free(negp);
		free(weight);
		free(pos_incov);
		free(neg_incov);

		for (i = 0; i < mpos; i++)
			free(poso[i]);

		for (i = 0; i < mneg; i++)
			free(nego[i]);

		free(poso);
		free(nego);

		is_train = false;
	}

	//Training
	void train(int DIM, int npos, int nneg, double **pos, double **neg);

	int nearestCent(double **cent, int nCent, int dim, double *pat);
	void getCentroid_LVQ(double **sample, double **result, double ***sigma, int n, int dim, int nCent, double learning_rate = 0.1, int epoch = 100);

	void clustering();
	void kernelizing();
	void lda();
	void beta();

	void getMisclassfiedSample();
	void getUncertainSample();

	bool getOutput(double* pat, double &val, double* &kval);

	double innerP(int dim, double *x, double *y) {
		int i;
		double val = 0.f;

		for (i = 0; i < dim; i++)
			val += x[i] * y[i];

		return val;
	}
	// All terms of cov. matrix
	double getKerOutCov(double *pat, int dim, double alpha, double *mean, double **incov);
	// Main diagonal terms of cov. matrix
	double getKerOut(double *pat, int dim, double alpha, double *mean, double **var);

	//instances
	bool is_train = false;
	bool is_libsvm;
	int DIM;
	int cpos, cneg;
	int npos, nneg;
	double alpha;
	double **pos, **neg;
	double **posm, **negm;
	double ***posv, ***negv;
	double ***pos_incov, ***neg_incov;
	double **posk, **negk;
	double *weight;
	double *posp, *negp;
	BD *posb, *negb;
	CPONBDA *bda;

	int mdim;
	int mpos, mneg;
	double **poso, **nego;
	double *maxo;
};
#endif