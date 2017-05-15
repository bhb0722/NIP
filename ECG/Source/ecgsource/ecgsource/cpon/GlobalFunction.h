/*
Copyright 2015. Bae Hanbin(bhb0722@nate.com or bhb0722@skku.edu) all rights reserved.
Beta Functions for Global Function
*/
#ifndef GLOBAL_FUN
#define GLOBAL_FUN

#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
#define LOG(tag, str) {cout << tag << " : " << str << endl;}
#define GET_PSUM for (j=0;j<ndim;j++){ for(sum=0.0,i=0;i<mpts;i++) sum+=p[i][j];\
				psum[j] = sum; }
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

#define NMAX 1000
#define TINY 1.0e-10
#define EPS 3.0e-7
#define FPMIN 1.0e-20
#define MAXIT 1000

double** getIterativeMatrixInversion(double **a, int n);
double** getMatrixInverse(double **a, int n);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);
// BetaCDF
// Beta function.
double getBetaCoeff(double x, double a, double b);
double getBetaCDF(double x, double params[]);

// Downhill Simplex Method
double getDWLS(double* sample, int n, double params[]);
void amoeba(double* sample, int n, double p[5][4], double y[], double ftol, int *nfunk);
double amotry(double* sample, int n, double p[5][4], double y[], double psum[], int ihi, double fac);

// FLDA
void getCentroid_LVQ(double **sample, double **result, double **sigma, int n, int dim, int nCent, double learning_rate = 0.1, int epoch = 100);
void getCentroid(double **sample, int n, double **result, double **sigma, int* K, int ndim);



#endif
