#pragma once
#ifndef CPON_1CENT_H
#define CPON_1CENT_H

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
class CPON_1CENT {


public:

	//Instructors
	CPON_1CENT() {
	}

	CPON_1CENT(double alpha) {
		this->alpha = alpha;
		TP = 0; TN = 0; FP = 0; FN = 0;
	}

	void clear() {
		int i;

		for (i = 0; i < npos; i++)
			free(posk[i]);

		for (i = 0; i < nneg; i++)
			free(negk[i]);

		for (i = 0; i < DIM; i++)
		{
			free(posv[i]);
			free(negv[i]);
			free(pos_incov[i]);
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
		free(pos_incov);
		free(neg_incov);
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
		getMisclassfiedSample();
		//getUncertainSample();
	}

	//Test
	void test_counting(int test_npos, int test_nneg, double **test_pos, double **test_neg) {
		int i, j, cpos = 1, cneg = 1, dim = 2;
		double pr, p, n;
		double **test_posk, **test_negk;
		double *test_posp, *test_negp;

		test_posk = (double**)malloc(sizeof(double*) * test_npos);
		test_negk = (double**)malloc(sizeof(double*) * test_nneg);

		for (i = 0; i < test_npos; i++) {
			test_posk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

			for (j = 0; j < cpos; j++)
				test_posk[i][j] = getKerOutCov(pos[i], DIM, alpha, posm, pos_incov);

			for (j = 0; j < cneg; j++)
				test_posk[i][j + cpos] = getKerOutCov(pos[i], DIM, alpha, negm, neg_incov);
		}

		for (i = 0; i < test_nneg; i++) {
			test_negk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

			for (j = 0; j < cpos; j++)
				test_negk[i][j] = getKerOutCov(neg[i], DIM, alpha, posm, pos_incov);

			for (j = 0; j < cneg; j++)
				test_negk[i][j + cpos] = getKerOutCov(neg[i], DIM, alpha, negm, neg_incov);
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
		fout << "Precision," << pre1 << endl;
		fout << "Recall," << rec1 << endl;
		fout << "F1," << f11 << endl;
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

	double *posm, *negm;
	double **posv, **negv;
	double **pos_incov, **neg_incov;

	double **posk, **negk;
	double *weight;
	double *posp, *negp;

	int mpos, mneg;
	double **poso, **nego;

	BD *posb, *negb;
	CPONBDA *bda;

	int TP, TN, FP, FN;
	double acc, pre1, pre2, rec1, rec2, f11, f12;

	bool getOutput(double* pat, double &val, double* &kval) {
		double wval;
		double p, n;
		
		kval[0] = getKerOutCov(pat, DIM, alpha, posm, pos_incov);
		kval[1] = getKerOutCov(pat, DIM, alpha, negm, neg_incov);
		
		wval = innerP(2, kval, weight);

		p = posb->getCDFValue(wval);
		n = 1.f - negb->getCDFValue(wval);
		val = p / (p + n);

		return bda->isUncertain(wval);
	}

	void getMisclassfiedSample() {
		int i, j, dim = 2; // # of cent = 2
		double val, *kval = (double*)malloc(sizeof(double)*2);
		vector<int> tlab;

		//			posp[i] = innerP(dim, weight, posk[i]);

		for (i = 0; i < npos; i++) {
			getOutput(pos[i], val, kval);
			if (val < 0.5f) {
				tlab.push_back(i);
			}
		}

		mpos = tlab.size();

		poso = (double**)malloc(sizeof(double*) * mpos);

		for (i = 0; i < mpos; i++) {
			poso[i] = (double*)malloc(sizeof(double) * dim);

			for (j = 0; j < dim; j++)
				poso[i][j] = posk[tlab[i]][j];
		}

		tlab.clear();
		tlab.resize(0);

		for (i = 0; i < nneg; i++) {
			getOutput(neg[i], val, kval);
			val = 1.f - val;
			if (val < 0.5f) {
				tlab.push_back(i);
			}
		}

		mneg = tlab.size();

		nego = (double**)malloc(sizeof(double*) * mneg);

		for (i = 0; i < mneg; i++) {
			nego[i] = (double*)malloc(sizeof(double) * dim);

			for (j = 0; j < dim; j++)
				nego[i][j] = negk[tlab[i]][j];
		}

#ifdef PRINT_HIST
		ofstream fout("normal_uncertain.txt");

		fout << mpos << ' ' << dim << endl;
		for (i = 0; i < mpos; i++) {
			for (j = 0; j < dim; j++)
				fout << poso[i][j] << ' ';
			fout << endl;
		}

		fout.close();

		fout.open("abnormal_uncertain.txt");
		fout << mneg << ' ' << dim << endl;
		for (i = 0; i < mneg; i++) {
			for (j = 0; j < dim; j++)
				fout << nego[i][j] << ' ';
			fout << endl;
		}
		fout.close();
#endif
	}

	void getUncertainSample() {
		int i, j, dim = 2; // # of cent = 2
		double val;
		vector<int> tlab;

		//			posp[i] = innerP(dim, weight, posk[i]);

		for (i = 0; i < npos; i++) {
			val = innerP(dim, weight, posk[i]);
			if (bda->isUncertain(val)) {
				tlab.push_back(i);
			}
		}

		mpos = tlab.size();

		poso = (double**)malloc(sizeof(double*) * mpos);

		for (i = 0; i < mpos; i++) {
			poso[i] = (double*)malloc(sizeof(double) * dim);

			for (j = 0; j < dim; j++)
				poso[i][j] = posk[tlab[i]][j];
		}

		tlab.clear();
		tlab.resize(0);

		for (i = 0; i < nneg; i++) {
			val = innerP(dim, weight, negk[i]);
			if (bda->isUncertain(val)) {
				tlab.push_back(i);
			}
		}

		mneg = tlab.size();

		nego = (double**)malloc(sizeof(double*) * mneg);

		for (i = 0; i < mneg; i++) {
			nego[i] = (double*)malloc(sizeof(double) * dim);

			for (j = 0; j < dim; j++)
				nego[i][j] = negk[tlab[i]][j];
		}

#ifdef PRINT_HIST
		ofstream fout("normal_uncertain.txt");

		fout << mpos << ' ' << dim << endl;
		for (i = 0; i < mpos; i++) {
			for (j = 0; j < dim; j++)
				fout << poso[i][j] << ' ';
			fout << endl;
		}

		fout.close();

		fout.open("abnormal_uncertain.txt");
		fout << mneg << ' ' << dim <<  endl;
		for (i = 0; i < mneg; i++) {
			for (j = 0; j < dim; j++)
				fout << nego[i][j] << ' ';
			fout << endl;
		}
		fout.close();
#endif

	}

	//protected: // Training functions
	double getKerOutCov(double *pat, int dim, double alpha, double *mean, double **incov) {
		int i, j;
		double sum = 0.f, flo = 1.e-10f;
		double *temp = (double*) malloc(sizeof(double) * dim);

		for (i = 0; i < dim; i++)
			temp[i] = 0.f;

		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++)
				temp[i] += (pat[i] - mean[i]) * incov[j][i];
		}

		for(i=0;i<dim;i++)
			sum += temp[i] * (pat[i] - mean[i]);

		free(temp);

		sum *= -0.5f * alpha;
		sum = exp(sum);

		if (sum < flo)
			return flo;

		return sum;
	}

	void clustering() {
		int i,j,k;
		
		posm = (double*)malloc(sizeof(double) * DIM);
		negm = (double*)malloc(sizeof(double) * DIM);
		posv = (double**)malloc(sizeof(double*) * DIM);
		negv = (double**)malloc(sizeof(double*) * DIM);

		for (i = 0; i < DIM; i++)
		{
			posv[i] = (double*)malloc(sizeof(double) * DIM);
			negv[i] = (double*)malloc(sizeof(double) * DIM);
		}

		for (i = 0; i < DIM; i++) {
			posm[i] = 0.f;
			negm[i] = 0.f;

			for (j = 0; j < DIM; j++) {
				posv[i][j] = 0.f;
				negv[i][j] = 0.f;
			}
		}

		for (i = 0; i < npos; i++) {
			for (j = 0; j < DIM; j++)
				posm[j] += pos[i][j];
		}

		for (i = 0; i < DIM; i++)
			posm[i] /= (double)npos;

		for (i = 0; i < nneg; i++) {
			for (j = 0; j < DIM; j++)
				negm[j] += neg[i][j];
		}

		for (i = 0; i < DIM; i++)
			negm[i] /= (double)nneg;

		for (i = 0; i < npos; i++) {
			for (j = 0; j < DIM; j++) {
				for (k = 0; k < DIM; k++)
					posv[j][k] += (pos[i][j] - posm[j]) * (pos[i][k] - posm[k]);
			}
		}

		for (i = 0; i < nneg; i++) {
			for (j = 0; j < DIM; j++) {
				for (k = 0; k < DIM; k++)
					negv[j][k] += (neg[i][j] - negm[j]) * (neg[i][k] - negm[k]);
			}
		}

		for (i = 0; i < DIM; i++) {
			for (j = 0; j < DIM; j++) {
				posv[i][j] /= (double)npos;
				negv[i][j] /= (double)nneg;
			}
		}

		pos_incov = getIterativeMatrixInversion(posv,DIM);
		neg_incov = getIterativeMatrixInversion(negv,DIM);
	}

	void kernelizing() {
		int i, j, cpos = 1, cneg = 1;

		posk = (double**)malloc(sizeof(double*) * npos);
		negk = (double**)malloc(sizeof(double*) * nneg);

		for (i = 0; i < npos; i++) {
			posk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

			for (j = 0; j < cpos; j++)
				posk[i][j] = getKerOutCov(pos[i], DIM, alpha, posm, pos_incov);

			for (j = 0; j < cneg; j++)
				posk[i][j + cpos] = getKerOutCov(pos[i], DIM, alpha, negm, neg_incov);
		}

		for (i = 0; i < nneg; i++) {
			negk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

			for (j = 0; j < cpos; j++)
			{
				if (j == 64)
					j = j;
				negk[i][j] = getKerOutCov(neg[i], DIM, alpha, posm, pos_incov);

			}
			for (j = 0; j < cneg; j++)
				negk[i][j + cpos] = getKerOutCov(neg[i], DIM, alpha, negm, neg_incov);
		}
	}

	double innerP(int dim, double *x, double *y) {
		int i;
		double val = 0.f;

		for (i = 0; i < dim; i++)
			val += x[i] * y[i];

		return val;
	}

	void lda() {
		int i;
		int dim = 2;
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
		bda = new CPONBDA(posp, negp, npos, nneg);

		posb = bda->pos;
		negb = bda->neg;

#ifdef PRINT_PARAM
		ofstream fout("normal_param.txt");
		fout << posb->mParam[0] << endl;
		fout << posb->mParam[1] << endl;
		fout << posb->mParam[2] << endl;
		fout << posb->mParam[3] << endl;
		fout.close();
		fout.open("abnormal_param.txt");
		fout << negb->mParam[0] << endl;
		fout << negb->mParam[1] << endl;
		fout << negb->mParam[2] << endl;
		fout << negb->mParam[3] << endl;
		fout.close();
#endif
	}


};
#endif