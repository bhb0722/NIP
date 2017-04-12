#ifndef KNN_H
#define KNN_H

#pragma once
#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;
class kNN {
public:


	//Instructors
	kNN() {

	}

	kNN(int k) {
		this->k = k;

		TP = 0; TN = 0; FP = 0; FN = 0;
	}

	void clear() {
		free(map);
	}

	//Training
	void train(int dim, int npos, int nneg, double **pos, double **neg) {
		int i;

		this->dim = dim;
		this->npos = npos;
		this->nneg = nneg;
		this->pos = pos;
		this->neg = neg;

		n = npos + nneg;
		map = (double**)malloc(sizeof(double) * n);

		for (i = 0; i < npos; i++)
			map[i] = pos[i];

		for (i = 0; i < nneg; i++)
			map[i + npos] = neg[i];
	}

	double sqDist(double *x, double *y) {
		int i;
		double sum = 0;

		for (i = 0; i < dim; i++)
			sum += (x[i] - y[i])*(x[i] - y[i]);

		return sum;
	}

	bool isPos(double *pat) {
		int i, j, l, pc = 0, nc;
		int *ind = (int*)malloc(sizeof(int) * k);
		double *dist = (double*)malloc(sizeof(double) * k), idist;

		dist[0] = -1.f;

		for (i = 0; i < n; i++) {
			idist = sqDist(pat, map[i]);

			for (j = 0; j < k; j++) {
				if (dist[j] == -1 || idist < dist[j])
					break;
			}

			if (j < k) {
				for (l = k - 1; l > j; l--) {
					dist[l] = dist[l - 1];
					ind[l] = ind[l - 1];
				}

				dist[j] = idist;
				ind[j] = i;
			}
		}

		for (i = 0; i < k; i++) {
			if (ind[i] < npos)
				pc++;
		}

		nc = k - pc;

		return pc > nc;
	}

	//Test
	void test_counting(int test_npos, int test_nneg, double **test_pos, double **test_neg) {
		int i;

		for (i = 0; i < test_npos; i++) {
			if (isPos(test_pos[i])) TP++;
			else FP++;
		}

		for (i = 0; i < test_nneg; i++) {
			if (isPos(test_neg[i])) FN++;
			else TN++;
		}
	}

	void test_measure(char *fname) {
		char fch[100];
		ofstream fout;

		sprintf(fch, "%s_knn.csv", fname);
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
		cout << "Precision (NR) : " << pre1 << endl;
		cout << "Recall : " << rec1 << endl;
		cout << "F1 : " << f11 << endl;
		cout << "---------------AN----------------" << endl;
		cout << "Precision (NR) : " << pre2 << endl;
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

	int n, dim;
	double **map;
	int k;
	int npos, nneg;
	double **pos, **neg;

	int TP, TN, FP, FN;
	double acc, pre1, pre2, rec1, rec2, f11, f12;
};

class kNN_CENT {
public:


	//Instructors
	kNN_CENT() {
	}

	kNN_CENT(int k, int nCentPos, int nCentNeg) {
		this->k = k;
		this->nCentPos = nCentPos;
		this->nCentNeg = nCentNeg;

		TP = 0; TN = 0; FP = 0; FN = 0;
	}

	void clear() {
		int i, cpos = nCentPos, cneg = nCentNeg;

		free(map);

		for (i = 0; i < cpos; i++) {
			free(posm[i]);
			free(posv[i]);
		}

		for (i = 0; i<cneg; i++) {
			free(negm[i]);
			free(negv[i]);
		}

		free(posm);
		free(posv);
		free(negm);
		free(negv);
	}

	//Training
	void train(int dim, int npos, int nneg, double **pos, double **neg) {
		int i;

		this->dim = dim;
		this->npos = npos;
		this->nneg = nneg;
		this->pos = pos;
		this->neg = neg;

		clustering();

		n = nCentPos + nCentNeg;
		map = (double**)malloc(sizeof(double) * n);
		
		for (i = 0; i < nCentPos; i++)
			map[i] = posm[i];

		for (i = 0; i < nCentNeg; i++)
			map[i + nCentPos] = negm[i];
	}

	double sqDist(double *x, double *y) {
		int i;
		double sum = 0;

		for (i = 0; i < dim; i++)
			sum += (x[i] - y[i])*(x[i] - y[i]);

		return sum;
	}

	bool isPos(double *pat) {
		int i, j, l, pc = 0, nc;
		int *ind = (int*)malloc(sizeof(int) * k);
		double *dist = (double*)malloc(sizeof(double) * k), idist;

		dist[0] = -1.f;

		for (i = 0; i < n; i++) {
			idist = sqDist(pat, map[i]);

			for (j = 0; j < k; j++) {
				if (dist[j] == -1 || idist < dist[j])
					break;
			}

			if (j < k) {
				for (l = k-1; l > j; l--) {
					dist[l] = dist[l - 1];
					ind[l] = ind[l - 1];
				}

				dist[j] = idist;
				ind[j] = i;
			}
		}

		for (i = 0; i < k; i++) {
			if (ind[i] < nCentPos) 
				pc++;
		}

		nc = k - pc;

		return pc > nc;
	}

	//Test
	void test_counting(int test_npos, int test_nneg, double **test_pos, double **test_neg) {
		int i;
		
		for (i = 0; i < test_npos; i++) {
			if (isPos(test_pos[i])) TP++;
			else FP++;
		}

		for (i = 0; i < test_nneg; i++) {
			if (isPos(test_neg[i])) FN++;
			else TN++;
		}
	}

	void test_measure(char *fname) {
		char fch[100];
		ofstream fout;

		sprintf(fch, "%s_knn.csv", fname);
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
		cout << "Precision (NR) : " << pre1 << endl;
		cout << "Recall : " << rec1 << endl;
		cout << "F1 : " << f11 << endl;
		cout << "---------------AN----------------" << endl;
		cout << "Precision (NR) : " << pre2 << endl;
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

	void clustering() {
		int i;
		int cpos = nCentPos, cneg = nCentNeg;
		// void getCentroid(double **sample, int n, double **result, double **sigma, int* K, int ndim)

		posm = (double**)malloc(sizeof(double*) * cpos);
		posv = (double**)malloc(sizeof(double*) * cpos);
		negm = (double**)malloc(sizeof(double*) * cneg);
		negv = (double**)malloc(sizeof(double*) * cneg);

		for (i = 0; i < cpos; i++) {
			posm[i] = (double*)malloc(sizeof(double) * dim);
			posv[i] = (double*)malloc(sizeof(double) * dim);
		}

		for (i = 0; i<cneg; i++) {
			negm[i] = (double*)malloc(sizeof(double) * dim);
			negv[i] = (double*)malloc(sizeof(double) * dim);
		}

		getCentroid(pos, npos, posm, posv, &cpos, dim);
		getCentroid(neg, nneg, negm, negv, &cneg, dim);
	}


	//instances

	int n, dim;
	double **map;
	int k;
	int npos, nneg;
	double **pos, **neg;
	int nCentPos, nCentNeg;
	double **posm, **posv, **negm, **negv;

	int TP, TN, FP, FN;
	double acc, pre1, pre2, rec1, rec2, f11, f12;
};

#endif