#include <iostream>
#include <fstream>
#include <cstdio>
#include "CPON.h"
#include "kNN.h"
#include "svm.h"
using namespace std;

enum feattype {
	alskew,
	alkur,
	alsd,
	alrm,
	sdnn,
	rmssd,
	sdrm,
	al,
	mvrrd,
	alvrrd,
};


// setting factor
#define TYPE al


// K-Fold Validation
#define KFOLD 10

// CPON Parameters
#define nCentPos 220
#define nCentNeg 110
#define ALPHA 0.5f
#define RDIM 9

// k-Nearest Neighbor Parameters
#define K 1

// 

int DIM;
int npos, nneg;
double **rawp, **rawn;
void input();
void feature(int type);
void kfold();


//	alskew, alkur, alsd, alrm, sdnn, rmssd, sdrm
//	loga, logl, skew, kurto, mean, SDNN, RMSSD;
double* setfeat(double *pat, int type) {
	double* feat = (double*)malloc(sizeof(double) * DIM);

	switch (type) {
	case alskew:
		//DIM = 3;
		feat[0] = pat[0];
		feat[1] = pat[1];
		feat[2] = pat[2];
		break;
	case alkur:
		//DIM = 3;
		feat[0] = pat[0];
		feat[1] = pat[1];
		feat[2] = pat[3];
		break;
	case alsd:
		//DIM = 3;
		feat[0] = pat[0];
		feat[1] = pat[1];
		feat[2] = pat[5];
		break;
	case alrm:
		//DIM = 3;
		feat[0] = pat[0];
		feat[1] = pat[1];
		feat[2] = pat[6];
		break;
	case sdnn:
		//DIM = 1;
		feat[0] = pat[5];
		break;
	case rmssd:
		//DIM = 1;
		feat[0] = pat[6];
		break;
	case sdrm:
		//DIM = 2;
		feat[0] = pat[5];
		feat[1] = pat[6];
		break;
	case al:
		feat[0] = pat[0];
		feat[1] = pat[1];
		break;
	case alvrrd:
		feat[0] = pat[0];
		feat[1] = pat[1];
		feat[2] = pat[8];
		break;
	case mvrrd:
		feat[0] = pat[7];
		feat[1] = pat[8];
		break;
	}

	return feat;
}

void setdim(int type) {
	switch (type) {
	case alskew:
		DIM = 3;
		cout << "log a, log l, skewness" << endl;
		break;
	case alkur:
		DIM = 3;
		cout << "log a, log l, kurtosis" << endl;
		break;
	case alsd:
		DIM = 3;
		cout << "log a, log l, SDNN" << endl;
		break;
	case alrm:
		DIM = 3;
		cout << "log a, log l, RMSSD" << endl;
		break;
	case sdnn:
		DIM = 1;
		cout << "SDNN" << endl;
		break;
	case rmssd:
		DIM = 1;
		cout << "RMSSD" << endl;
		break;
	case sdrm:
		DIM = 2;
		cout << "SDNN, RMSSD" << endl;
		break;
	case al:
		DIM = 2;
		cout << "log a, log l" << endl;
		break;
	case alvrrd:
		DIM = 3;
		cout << "log a, log l, var. of diff. of RR int.s" << endl;
		break;
	case mvrrd:
		DIM = 2;
		cout << "Mean, Var. of diff. of RR int.s" << endl;
	}
}

int main() {
	input();
	setdim(TYPE);
	kfold();
	return 0;
}


void input() {
	int n, i, j, lab;
	ifstream fin("normal.txt");

	fin >> n;
	npos = n;

	rawp = (double**)malloc(sizeof(double*) * n);

	for (i = 0; i < n; i++) {
		rawp[i] = (double*)malloc(sizeof(double) * RDIM);

		for (j = 0; j < RDIM; j++)
			fin >> rawp[i][j];
	}

	fin.close();
	fin.open("abnormal.txt");

	fin >> n;
	nneg = n;
	
	rawn = (double**)malloc(sizeof(double) * n);

	for (i = 0; i < n; i++) {
		rawn[i] = (double*)malloc(sizeof(double) * RDIM);

		for (j = 0; j < RDIM; j++)
			fin >> rawn[i][j];

		fin >> lab;
	}

	fin.close();
}


void kfold() {
	int i, tt_j, tr_j, k, st, ed;
	int tr_npos, tr_nneg;
	double **tr_pos, **tr_neg;
	int tt_npos, tt_nneg;
	double **tt_pos, **tt_neg;
	double factor = 1.f / KFOLD;
	CPON *cpon = new CPON(nCentPos, nCentNeg, ALPHA);
	kNN *knn = new kNN(1);
	SVM *svm = new SVM();
	//kNN_CENT *knn = new kNN_CENT(K, nCentPos, nCentNeg);

	for (k = 0; k < KFOLD; k++) {
		//////////////////////POS
		st = k * factor * npos;
		ed = (k + 1) * factor * npos;
		tt_npos = ed - st;
		tr_npos = npos - tt_npos;

		tt_pos = (double**)malloc(sizeof(double*) * tt_npos);
		tr_pos = (double**)malloc(sizeof(double*) * tr_npos);
		tt_j = 0;
		tr_j = 0;

		for (i = 0; i < npos; i++) {
			if (i >= st && i < ed) //test
				tt_pos[tt_j++] = setfeat(rawp[i], TYPE);
			else
				tr_pos[tr_j++] = setfeat(rawp[i], TYPE);
		}

		//////////////////////NEG
		st = k * factor * nneg;
		ed = (k + 1) * factor * nneg;
		tt_nneg = ed - st;
		tr_nneg = nneg - tt_nneg;

		tt_neg = (double**)malloc(sizeof(double*) * tt_nneg);
		tr_neg = (double**)malloc(sizeof(double*) * tr_nneg);
		tt_j = 0;
		tr_j = 0;

		for (i = 0; i < nneg; i++) {
			if (i >= st && i < ed) //test
				tt_neg[tt_j++] = setfeat(rawn[i], TYPE);
			else
				tr_neg[tr_j++] = setfeat(rawn[i], TYPE);
		}


		cpon->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		cpon->test_counting(tt_npos, tt_nneg, tt_pos, tt_neg);
		cpon->clear();

		knn->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		knn->test_counting(tt_npos, tt_nneg, tt_pos, tt_neg);
		knn->clear();
		
		svm->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		svm->test_counting(tt_npos, tt_nneg, tt_pos, tt_neg);
		svm->clear();

		/////////////////////CLEAR
		for (i = 0; i < tr_npos; i++) free(tr_pos[i]);
		for (i = 0; i < tr_nneg; i++) free(tr_neg[i]);
		for (i = 0; i < tt_npos; i++) free(tt_pos[i]);
		for (i = 0; i < tt_nneg; i++) free(tt_neg[i]);
		free(tr_pos);
		free(tr_neg);
		free(tt_pos);
		free(tt_neg);
	}

	cpon->test_measure("result");
	knn->test_measure("result");
	svm->test_measure("result");
}



