#ifndef CPON_LIBSVM_H
#define CPON_LIBSVM_H

#ifndef CPON_H
#include "cpon/betadist.h"
#include "cpon/GlobalFunction.h"
#endif
#include "svm.h"
#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;

#define PRINT_CENT 0 // 금방됨
#define PRINT_KERNEL 0 // 금방됨
#define PRINT_HIST 0 // 조금걸림
#define PRINT_PARAM 0 // 걸림
#pragma warning(disable:4996)
class CPON_LIBSVM {


public:

	//Instructors
	CPON_LIBSVM() {
		TP = 0; TN = 0; FP = 0; FN = 0;
	}

	void clear() {
		
	}

	/*
	@param
		gamma : kernel width
	@ret
		pvalue of positive : p_1
		pvalue of negative : p_2
		return p_1 * p_2
	*/
	double train_libsvm(SVM *svm) {
		int i, j;
		double p_1, p_2;
		svm_problem *prob;
		svm_node *inode;
		svm_model *model;

		//svm->train(DIM, npos, nneg, pos, neg);
		prob = &(svm->prob);
		model = svm->model;

		posp = (double*)malloc(sizeof(double) * npos);
		negp = (double*)malloc(sizeof(double) * nneg);

		for(j = 0; j < npos; j++) {
			inode = prob->x[j];
			svm_predict_values(model, inode, &posp[j]);
		}

		for (i = 0; i < nneg; i++) {
			inode = prob->x[j++];
			svm_predict_values(model, inode, &negp[i]);
		}

		posb = new BD(&posp[0], npos);
		negb = new BD(&negp[0], nneg);

		posb->getEstimatorbyCG();
		negb->getEstimatorbyCG();
		//posb->getEstimator();
		//negb->getEstimator();
		p_1 = posb->pValue;
		p_2 = negb->pValue;

		free(posp);
		free(negp);

		return  p_1 * p_2;
	}

	//Training
	void train(SVM *svm, double gamma, int npos, int nneg) {
		int j;
		double pv;

		this->npos = npos;
		this->nneg = nneg;
		
		posb = NULL;
		negb = NULL;

		pv = train_libsvm(svm);
			
		if (max_pv < pv){
			for (j = 0; j < 4; j++) {
				par1[j] = posb->mParam[j];
				par2[j] = negb->mParam[j];
			}
			
			max_pv = pv;
			max_g = gamma;
		}
	}


	//Test
	void test_counting(SVM *svm, int test_npos, int test_nneg, double **test_pos, double **test_neg) {
		int i;
		double o, pr, p, n;

		posb = new BD(npos, npos, par1);
		negb = new BD(nneg,nneg,par2);

		for (i = 0; i < test_npos; i++) { // if pr >= 0.5f then TP, otherwise FP
										  //getKerOut(pos[i], DIM, alpha, posm[j], posv[j]);
			o = svm->getOutput(test_pos[i]);
			p = posb->getCDFValue(o);
			n = 1.f - negb->getCDFValue(o);
			pr = p / (p + n);

			if (pr >= 0.5f) TP++;
			else FP++;
		}

		for (i = 0; i < test_nneg; i++) { // if pr >= 0.5f then TN, otherwise FN
			o = svm->getOutput(test_neg[i]);
			p = posb->getCDFValue(o);
			n = 1.f - negb->getCDFValue(o);
			pr = n / (p + n);

			if (pr >= 0.5f) TN++;
			else FN++;
		}
	}

	void test_measure(char *fname) {
		char fch[100];
		ofstream fout;

		sprintf(fch, "%s_cpon_svm.csv", fname);
		fout.open(fch);
		acc = 1.f - ((double)(FP + FN)) / (TP + FP + TN + FN);
		pre1 = (double)TP / (TP + FP);
		rec1 = (double)TP / (TP + FN);
		f11 = (2.f * pre1 * rec1) / (pre1 + rec1);
		pre2 = (double)TN / (TN + FN);
		rec2 = (double)TN / (TN + FP);
		f12 = (2.f * pre2 * rec2) / (pre2 + rec2);

		cout << "--------------CPON---------------" << endl;
		cout << "Accuracy : " << acc << endl;
		cout << "---------------NR----------------" << endl;
		cout << "Precision : " << pre1 << endl;
		cout << "Recall : " << rec1 << endl;
		cout << "F1 : " << f11 << endl;
		cout << "---------------AN----------------" << endl;
		cout << "Precision : " << pre2 << endl;
		cout << "Recall : " << rec2 << endl;
		cout << "F1 : " << f12 << endl;
		cout << endl;


		fout << "Accuracy," << acc << endl;
		fout << "Precision," << pre1 << endl;
		fout << "Recall," << rec1 << endl;
		fout << "F1," << f11 << endl;
		fout << "Precision," << pre2 << endl;
		fout << "Recall," << rec2 << endl;
		fout << "F1," << f12 << endl;
		fout.close();
	}

	//instances

	int DIM;
	int npos, nneg;
	double alpha;
	double **pos, **neg;
	double *posp, *negp;
	double max_pv = 0.f;
	double max_g;
	double par1[4], par2[4];
	BD *posb, *negb;
	int TP, TN, FP, FN;
	double acc, pre1, pre2, rec1, rec2, f11, f12;
};
#endif