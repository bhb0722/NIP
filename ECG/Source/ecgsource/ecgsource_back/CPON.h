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
#define PRINT_HIST 1 // 조금걸림
#define PRINT_PARAM 0 // 걸림
#pragma warning(disable:4996)
class CPON {


public:

	//Instructors
	CPON() {
	}

	CPON(int nCentPos, int nCentNeg, double alpha) {
		this->nCentPos = nCentPos;
		this->nCentNeg = nCentNeg;
		this->alpha = alpha;
		TP = 0; TN = 0; FP = 0; FN = 0;
	}

	void clear() {
		int i;

		for (i = 0; i < npos; i++)
			free(posk[i]);

		for (i = 0; i < nneg; i++)
			free(negk[i]);

		for (i = 0; i < nCentPos; i++)
		{
			free(posm[i]);
			free(posv[i]);
		}

		for (i = 0; i < nCentNeg; i++) {
			free(negm[i]);
			free(negv[i]);
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
	}

	//Training
	void train(int DIM, int npos, int nneg, double **pos, double **neg) {
		this->DIM = DIM;
		this->npos = npos;
		this->nneg = nneg;
		this->pos = pos;
		this->neg = neg;

		clustering();
		kernelizing();
		lda();
		beta();
	}

	//Test
	void test_counting(int test_npos, int test_nneg, double **test_pos, double **test_neg) {
		int i, j, cpos = nCentPos, cneg = nCentNeg, dim = nCentPos + nCentNeg;
		double pr, p, n;
		double **test_posk, **test_negk;
		double *test_posp, *test_negp;

		test_posk = (double**)malloc(sizeof(double*) * test_npos);
		test_negk = (double**)malloc(sizeof(double*) * test_nneg);

		for (i = 0; i < test_npos; i++) {
			test_posk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

			for (j = 0; j < cpos; j++)
				test_posk[i][j] = getKerOut(test_pos[i], DIM, alpha, posm[j], posv[j]); //posm , posv

			for (j = 0; j < cneg; j++)
				test_posk[i][j + cpos] = getKerOut(test_pos[i], DIM, alpha, negm[j], negv[j]); //negm, negv
		}

		for (i = 0; i < test_nneg; i++) {
			test_negk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

			for (j = 0; j < cpos; j++)
				test_negk[i][j] = getKerOut(test_neg[i], DIM, alpha, posm[j], posv[j]); //posm, posv

			for (j = 0; j < cneg; j++)
				test_negk[i][j + cpos] = getKerOut(test_neg[i], DIM, alpha, negm[j], negv[j]); //negm, negv
		}

		test_posp = (double*)malloc(sizeof(double) * test_npos);
		test_negp = (double*)malloc(sizeof(double) * test_nneg);

		for (i = 0; i < test_npos; i++)
			test_posp[i] = innerP(dim, weight, test_posk[i]);

		for (i = 0; i < test_nneg; i++)
			test_negp[i] = innerP(dim, weight, test_negk[i]);

		for (i = 0; i < test_npos; i++) { // if pr >= 0.5f then TP, otherwise FP
										  //getKerOut(pos[i], DIM, alpha, posm[j], posv[j]);
			p = posb->getCDFValue(test_posp[i]);
			n = 1.f - negb->getCDFValue(test_posp[i]);
			pr = p / (p + n);

			if (pr >= 0.5f) TP++;
			else FP++;

			free(test_posk[i]);
		}

		for (i = 0; i < test_nneg; i++) { // if pr >= 0.5f then TN, otherwise FN
			p = posb->getCDFValue(test_negp[i]);
			n = 1.f - negb->getCDFValue(test_negp[i]);
			pr = n / (p + n);

			if (pr >= 0.5f) TN++;
			else FN++;

			free(test_negk[i]);
		}

		free(test_posp);
		free(test_negp);
	}

	void test_measure(char *fname) {
		char fch[100];
		ofstream fout;

		sprintf(fch, "%s_cpon.csv", fname);
		fout.open(fch);
		acc = 1.f - ((double)(FP + FN)) / (TP + FP + TN + FN);
		pre1 = (double)TP / (TP + FP);
		rec1 = (double)TP / (TP + FN);
		f11 = (2.f * pre1 * rec1) / (pre1 + rec1);
		pre2 = (double)TN / (TN + FN);
		rec2 = (double)TN / (TN + FP);
		f12 = (2.f * pre2 * rec2) / (pre2 + rec2);


		cout << "Accuracy : " << acc << endl;
		cout << "---------------NR----------------" << endl;
		cout << "Precision : " << pre1 << endl;
		cout << "Recall : " << rec1 << endl;
		cout << "F1 : " << f11 << endl;
		cout << "---------------AN----------------" << endl;
		cout << "Precision : " << pre2 << endl;
		cout << "Recall : " << rec2 << endl;
		cout << "F1 : " << f12 << endl;

		
		fout << "Accuracy," << acc << endl;
		fout << "Normal" << endl;
		fout << "Precision," << pre1 << endl;
		fout << "Recall," << rec1 << endl;
		fout << "F1," << f11 << endl;
		fout << "Abnormal" << endl;
		fout << "Precision," << pre2 << endl;
		fout << "Recall," << rec2 << endl;
		fout << "F1," << f12 << endl; 
	}

	//instances

	int DIM;
	int nCentPos, nCentNeg;
	int npos, nneg;
	double alpha;
	double **pos, **neg;
	double **posm, **negm;
	double **posv, **negv;
	double **posk, **negk;
	double *weight;
	double *posp, *negp;
	BD *posb, *negb;
	int TP, TN, FP, FN ;
	double acc, pre1, pre2, rec1, rec2, f11, f12;

//protected: // Training functions
	double getKerOut(double *pat, int dim, double alpha, double *mean, double *var) {
		int i;
		double sum = 0.f, flo = 1.e-20f;

		for (i = 0; i < dim; i++)
			sum += (pat[i] - mean[i]) * (pat[i] - mean[i]) / var[i];

		sum *= -0.5f * alpha;
		sum = exp(sum);

		if (sum < flo)
			return flo;

		return sum;
	}

	void clustering() {
		int i;
		int cpos = nCentPos, cneg = nCentNeg;
		// void getCentroid(double **sample, int n, double **result, double **sigma, int* K, int ndim)

		posm = (double**)malloc(sizeof(double*) * cpos);
		posv = (double**)malloc(sizeof(double*) * cpos);
		negm = (double**)malloc(sizeof(double*) * cneg);
		negv = (double**)malloc(sizeof(double*) * cneg);

		for (i = 0; i < cpos; i++) {
			posm[i] = (double*)malloc(sizeof(double) * DIM);
			posv[i] = (double*)malloc(sizeof(double) * DIM);
		}

		for (i = 0; i<cneg; i++) {
			negm[i] = (double*)malloc(sizeof(double) * DIM);
			negv[i] = (double*)malloc(sizeof(double) * DIM);
		}

		getCentroid(pos, npos, posm, posv, &cpos, DIM);
		getCentroid(neg, nneg, negm, negv, &cneg, DIM);
	}

	void kernelizing() {
		int i, j, cpos = nCentPos, cneg = nCentNeg;

		posk = (double**)malloc(sizeof(double*) * npos);
		negk = (double**)malloc(sizeof(double*) * nneg);

		for (i = 0; i < npos; i++) {
			posk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

			for (j = 0; j < cpos; j++)
				posk[i][j] = getKerOut(pos[i], DIM, alpha, posm[j], posv[j]); //posm , posv

			for (j = 0; j < cneg; j++)
				posk[i][j + cpos] = getKerOut(pos[i], DIM, alpha, negm[j], negv[j]); //negm, negv
		}

		for (i = 0; i < nneg; i++) {
			negk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

			for (j = 0; j < cpos; j++)
				negk[i][j] = getKerOut(neg[i], DIM, alpha, posm[j], posv[j]); //posm, posv

			for (j = 0; j < cneg; j++)
				negk[i][j + cpos] = getKerOut(neg[i], DIM, alpha, negm[j], negv[j]); //negm, negv
		}
	}

	double innerP(int dim, double *x, double *y) {
		int i;
		double val = 0.f;

		for (i = 0; i < dim; i++)
			val += x[i] * y[i];

		return -1.f * val;
	}

	void lda() {
		int i;
		int dim = nCentPos + nCentNeg;
		FLDA flda;

		weight = (double*)malloc(sizeof(double) * dim);

		flda.weight = weight;
		flda.nPos = npos;
		flda.nNeg = nneg;
		flda.getFisherCriterion(posk, negk, dim);

		posp = (double*)malloc(sizeof(double) * npos);
		negp = (double*)malloc(sizeof(double) * nneg);

		for (i = 0; i < npos; i++)
			posp[i] = innerP(dim, weight, posk[i]);

		for (i = 0; i < nneg; i++)
			negp[i] = innerP(dim, weight, negk[i]);

#ifdef PRINT_HIST
		ofstream fout("normal_hist.txt");

		fout << npos << endl;
		for (i = 0; i < npos; i++)
			fout << posp[i] << endl;

		fout.close();

		fout.open("abnormal_hist.txt");
		fout << nneg << endl;
		for (i = 0; i < nneg; i++)
			fout << negp[i] << endl;
		fout.close();
#endif

		/*	for (int i = 0; i < dim; i++)
		cout << weight[i] << ' ';
		cout << endl;*/
	}

	void beta() {
		BD *pos, *neg;
		pos = new BD(posp, npos);
		neg = new BD(negp, nneg);
		pos->getEstimator();
		neg->getEstimator();
		posb = pos;
		negb = neg;

#ifdef PRINT_PARAM
		ofstream fout("normal_param.txt");
		fout << pos->mParam[0] << endl;
		fout << pos->mParam[1] << endl;
		fout << pos->mParam[2] << endl;
		fout << pos->mParam[3] << endl;
		fout.close();
		fout.open("abnormal_param.txt");
		fout << neg->mParam[0] << endl;
		fout << neg->mParam[1] << endl;
		fout << neg->mParam[2] << endl;
		fout << neg->mParam[3] << endl;
		fout.close();
#endif
	}


};
#endif