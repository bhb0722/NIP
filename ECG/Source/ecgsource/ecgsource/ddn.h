#pragma once
#ifndef DDN_H
#define DDN_H
#include "cpon.h"
#include "cpon_1cent.h"
#include <vector>
using namespace std;

class DDN {

public:
	int nlayer;
	double *alpha;

	int npos, nneg;
	double **pos, **neg;

	vector<CPON*> layers;

	int TP, TN, FP, FN;
	double acc, pre1, pre2, rec1, rec2, f11, f12;

	DDN(int nlayer, int* cpos, int* cneg, double* alpha) {
		int i;

		TP = TN = FP = FN = 0;

		this->nlayer = nlayer;
		this->alpha = alpha;

		layers.resize(nlayer);
		

		for (i = 0; i < nlayer; i++)
			layers[i] = new CPON(cpos[i],cneg[i], alpha[i]);
	}

	void train(int dim, int npos, int nneg, double **pos, double **neg) {
		int i;


		layers[0]->train(dim, npos, nneg, pos, neg);

		for (i = 1; i < nlayer; i++) {
			if(layers[i-1]->mpos >= 10 && layers[i-1]->mneg >= 10)
				layers[i]->train(layers[i - 1]->mdim, layers[i - 1]->mpos, layers[i - 1]->mneg, layers[i - 1]->poso, layers[i - 1]->nego);
		}
	}

	double getOutput(int dim, double* pat) {
		int i;
		double lout;
		double *in = pat, *out = NULL;

		for (i = 0; i < nlayer; i++) {
			if (i > 0) {
				if (i != 1)
					free(in);
				in = out;
			}

			out = (double*)malloc(sizeof(double) * layers[i]->mdim);
			if (layers[i]->is_train && !layers[i]->getOutput(in, lout, out))
				break;
		}

		free(out);

		return lout;
	}

	//Test
	void test_counting(int dim, int test_npos, int test_nneg, double **test_pos, double **test_neg) {
		int i;
		double pr;

		for (i = 0; i < test_npos; i++) { // if pr >= 0.5f then TP, otherwise FP
			pr = getOutput(dim, test_pos[i]);

			if (pr >= 0.5f) TP++;
			else FP++;
		}

		for (i = 0; i < test_nneg; i++) { // if pr >= 0.5f then TN, otherwise FN
			pr = 1.f - getOutput(dim, test_neg[i]);

			if (pr >= 0.5f) TN++;
			else FN++;
		}
	}

	void test_measure(char *fname) {
		char fch[100];
		ofstream fout;

		sprintf(fch, "%s_ddn.csv", fname);
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

	DDN() {

	}

	void clear() {
		int i;

		for (i = 0; i < nlayer; i++)
			layers[i]->clear();
	}
};

class DDN_1cent {

public:
	int nlayer;
	double *alpha;

	int npos, nneg;
	double **pos, **neg;

	vector<CPON_1CENT*> layers;

	int TP, TN, FP, FN;
	double acc, pre1, pre2, rec1, rec2, f11, f12;

	DDN_1cent(int nlayer,double* alpha) {
		int i;

		TP = TN = FP = FN = 0;

		this->nlayer = nlayer;
		this->alpha = alpha;

		layers.resize(nlayer);
		
		for (i = 0; i < nlayer; i++)
			layers[i] = new CPON_1CENT(alpha[i]);
	}

	void train(int dim, int npos, int nneg, double **pos, double **neg) {
		int i;

		layers[0]->train(dim, npos, nneg, pos, neg);
		
		for (i = 1; i < nlayer; i++)
			layers[i]->train(2, layers[i - 1]->mpos, layers[i - 1]->mneg, layers[i - 1]->poso, layers[i - 1]->nego);
	}

	double getOutput(int dim, double* pat) {
		int i;
		double lout;
		double *in = pat, *out = NULL;

		for (i = 0; i < nlayer; i++) {
			if (i > 0) {
				if(i != 1)
					free(in);
				in = out;
			}

			out = (double*)malloc(sizeof(double) * 2);
			if (!layers[i]->getOutput(in, lout, out))
				break;
		}

		free(out);

		return lout;
	}

	//Test
	void test_counting(int dim, int test_npos, int test_nneg, double **test_pos, double **test_neg) {
		int i;
		double pr;
		
		for (i = 0; i < test_npos; i++) { // if pr >= 0.5f then TP, otherwise FP
			pr = getOutput(dim, test_pos[i]);

			if (pr >= 0.5f) TP++;
			else FP++;
		}

		for (i = 0; i < test_nneg; i++) { // if pr >= 0.5f then TN, otherwise FN
			pr = 1.f - getOutput(dim, test_neg[i]);

			if (pr >= 0.5f) TN++;
			else FN++;
		}
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

	DDN_1cent() {

	}

	void clear() {
		int i;

		for (i = 0; i < nlayer; i++)
			layers[i]->clear();
	}
};

#endif