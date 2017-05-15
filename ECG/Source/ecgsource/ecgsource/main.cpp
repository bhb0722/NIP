#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include "CPON.h"
#include "CPON_libsvm.h"
#include "ddn.h"
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
int TYPE = alsd;
//#define TYPE al
#define BTN_CPON 1
#define BTN_KNN 0
#define BTN_SVM 0
#define BTN_CPON_SVM 0
#define BTN_CPON_1CENT 0
#define BTN_DDN 0


// K-Fold Validation
#define KFOLD 10

// CPON Parameters
#define nCentPos 22
#define nCentNeg 11
#define ALPHA 0.3
#define RDIM 9

// k-Nearest Neighbor Parameters
#define K 1

// CPON with SVM

int nG = 6;
double G[9] = {4.f, 2.f, 1.f, 0.5f, 0.25f, 0.125f};
int DIM;
int npos, nneg;
double **rawp, **rawn;
void input();
void feature(int type);
void kfold(char *fname);


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
	char fname[100];
	//srand(time(NULL));
	input();
	

	//for (int i = 0; i < 10; i++)
	//{
	//	TYPE = i;
		setdim(TYPE);
		sprintf(fname, "%d_", TYPE);
		kfold(fname);
	//}
	return 0;
}


void input() {
	int n, i, j, lab,ri;
	bool *tag;
	ifstream fin("normal.txt");

	fin >> n;
	npos = n;

	rawp = (double**)malloc(sizeof(double*) * n);
	tag = (bool*)calloc(n,sizeof(bool));

	for (i = 0; i < n; i++) {
		do {
			ri = ((double)rand() / (RAND_MAX + 1)) * n;
		} while (tag[ri]);
		tag[ri] = true;
		rawp[ri] = (double*)malloc(sizeof(double) * RDIM);

		for (j = 0; j < RDIM; j++)
			fin >> rawp[ri][j];
	}

	free(tag);
	fin.close();
	fin.open("abnormal.txt");

	fin >> n;
	nneg = n;
	tag = (bool*)calloc(n, sizeof(bool));

	
	rawn = (double**)malloc(sizeof(double) * n);

	for (i = 0; i < n; i++) {
		do {
			ri = ((double)rand() / (RAND_MAX + 1)) * n;
		} while (tag[ri]);
		tag[ri] = true;
		rawn[ri] = (double*)malloc(sizeof(double) * RDIM);

		for (j = 0; j < RDIM; j++)
			fin >> rawn[ri][j];

		fin >> lab;
	}

	free(tag);
	fin.close();
}


void kfold(char *fname) {
	int i, tt_j, tr_j, k, st, ed;
	int tr_npos, tr_nneg;
	double **tr_pos, **tr_neg;
	int tt_npos, tt_nneg;
	double **tt_pos, **tt_neg;
	double factor = 1.f / KFOLD;

#if BTN_CPON
	CPON *cpon = new CPON(nCentPos, nCentNeg, ALPHA);
#endif

#if BTN_KNN
	kNN *knn = new kNN(1);
#endif

#if BTN_SVM
	SVM *svm = new SVM();
#endif

#if BTN_CPON_SVM
	SVM *svm;
	CPON_LIBSVM *cpon_svm = new CPON_LIBSVM();
#endif

#if BTN_CPON_1CENT
	/*int nlayer = 1;
	double alphas[] = {1.1f};
	DDN_1cent *ddn = new DDN_1cent(nlayer, alphas);*/
	CPON_1CENT *cpon_1 = new CPON_1CENT(ALPHA);
	//CPON_1CENT *cpon_2 = new CPON_1CENT(ALPHA);
#endif 

	//kNN_CENT *knn = new kNN_CENT(K, nCentPos, nCentNeg);

#if BTN_DDN
	int nlayer = 2;
	int cpos[] = {1, 1};
	int cneg[] = {1, 1};
	double alphas[] = {0.9f, 1.2f};
	DDN *ddn = new DDN(nlayer, cpos, cneg, alphas);
#endif 

	for (k = 0; k < KFOLD; k++) {
		cout << "Doing " << k << " fold" << endl;
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

#if BTN_CPON
		cpon->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		cpon->test_counting(tt_npos, tt_nneg, tt_pos, tt_neg);
		cpon->clear();
#endif

#if BTN_KNN
		knn->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		knn->test_counting(tt_npos, tt_nneg, tt_pos, tt_neg);
		knn->clear();
#endif

#if BTN_SVM
		svm->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		svm->test_counting(tt_npos, tt_nneg, tt_pos, tt_neg);
		svm->clear();
#endif

#if BTN_CPON_SVM
		for (i = 0; i < nG; i++) {
			svm = new SVM(G[i]);
			svm->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
			cpon_svm->train(svm, G[i], tr_npos, tr_nneg);
			svm->clear();
			delete svm;
		}

		svm = new SVM(cpon_svm->max_g);
		svm->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		svm->test_counting(tt_npos, tt_nneg, tt_pos, tt_neg);
		cpon_svm->test_counting(svm, tt_npos, tt_nneg, tt_pos, tt_neg);

#endif

#if BTN_CPON_1CENT
		/*ddn->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		ddn->test_counting(2, tt_npos, tt_nneg, tt_pos, tt_neg);
		ddn->clear();*/

		cpon_1->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		cpon_1->test_counting(tt_npos, tt_nneg, tt_pos, tt_neg);
		cpon_1->clear();
#endif

#if BTN_DDN
		ddn->train(DIM, tr_npos, tr_nneg, tr_pos, tr_neg);
		ddn->test_counting(2, tt_npos, tt_nneg, tt_pos, tt_neg);
		ddn->clear();
#endif

		/////////////////////CLEAR
		for (i = 0; i < tr_npos; i++) 
			free(tr_pos[i]);
		for (i = 0; i < tr_nneg; i++) free(tr_neg[i]);
		for (i = 0; i < tt_npos; i++) free(tt_pos[i]);
		for (i = 0; i < tt_nneg; i++) free(tt_neg[i]);
		free(tr_pos);
		free(tr_neg);
		free(tt_pos);
		free(tt_neg);
	}

#if BTN_CPON
	cpon->test_measure("result");
#endif
#if BTN_KNN
	knn->test_measure("result");
#endif
#if BTN_SVM
	svm->test_measure("result");
#endif

#if BTN_CPON_SVM
	svm->test_measure(fname);
	cpon_svm->test_measure(fname);
	svm->clear();
#endif

#if BTN_CPON_1CENT
	//ddn->test_measure("result");
	cpon_1->test_measure("result");
#endif
	
#if BTN_DDN

	ddn->test_measure("result");
#endif
}



