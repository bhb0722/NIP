/*
	Copyright 2015. Bae Hanbin(bhb0722@nate.com or bhb0722@skku.edu) all rights reserved.
	Decision Network for ASR.
	Join The Project with ETRI(www.etri.re.kr).
*/
#ifndef BETADIST_H
#define BETADIST_H
#include "GlobalFunction.h"

#include <stdlib.h>
#include <math.h>


#include <sstream>
#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <algorithm>
using namespace std;
#define CONFLEVEL 0.01
#define SIMPSONSTEP 9999
#define CENTMIN 1
#define CENTMAX 20
#define UV 1

#define LOG(tag, str) {cout << tag << " : " << str << endl;}
#define ISKSTEST 1
#define ISINPUT 0
#define ISRELATED 0
#define ISFLDA 1
#define ISOUTPUT 

#define absPath "output/"
#define ext ".txt"
#define posStr "pos_m"
#define negStr "neg_m"
class CPONLAYER;
class CPONNODE;
class CPONBDA;
class BD;
class FLDA;

class CPONLAYER{
public:
	CPONLAYER(){};
	CPONLAYER(int sampleSize, int nodeSize, int *label, int *nPos, double **sample);
	CPONLAYER(int nodeSize, const char outDir[]);

	double* getClassOutput(double* input);
	double* getCPONOutput(double* input);
	//one shot
	void doTraining();

	void trainDNN();
	void trainDDN();

	//3 parts
	void trainInputClass();
	void trainInputClass2();
	void trainFLDA();
	void trainOutputClass();

	void saveData();
	void saveInput();
	void saveRelated();
	void saveFLDA();
	void saveOutput();
	double getOutput();
	
	CPONNODE **mNode;
	int nodeSize, sampleSize;
	int *label, *nPos;
	bool *bCheck;
	string outDir;

	double **sample;
	double **posSamples, **negSamples;

	vector<double> *vPos, *vNeg;
};

class CPONNODE{
public:
	CPONNODE(int iLabel);
	CPONNODE(int iLabel, int nodeSize, CPONBDA *input);
	CPONNODE(int iLabel, int nodeSize, double *posSample, double *negSample, int nPos, int nNeg);
	~CPONNODE();
	void setClassProb(double prob);
	void pushRelatedNode(int p);
	bool checkUncertainty(double x);
	double getOutput(double *pattern);
	void getOptFLDA(int ndim);
	void setTransDist();
	vector<int>* getDep();
	void pushUPos(double *pattern);
	void pushUNeg(double *pattern);
	void deleteUSamples();
	void sortRelatedNode();
	double getCPONOutput(double x, CPONBDA *ibda);
	double getOriginalCPONOutput(double x, CPONBDA *ibda);
	double getLogCPONOutput(double x, CPONBDA *ibda);
	double getDNNOutputCPON(double x);
	double getDNNOutputWeightedCPON(double x);
	double getDNNOutputExp2(double x);
	void getTransDist();


	FLDA **mFLDA;
	CPONBDA *mInput, *mOutput;
	int iLabel, ndim, opti;
	vector<int> *dep;
	
	bool bCertain, *bDep;
	double classProb;

	vector<double*> uPosSample, uNegSample;
	double *uTransPos, *uTransNeg;
};


class FLDA{
public:
	FLDA(){};
	FLDA(double **posSample, double **negSample, int nPos, int nNeg, int ndim, int nCent);
	FLDA(int ndim, int nCent, double **centPos, double **centNeg, double *sigPos, double *sigNeg, double *weight);
	virtual ~FLDA();
	
	double getJ();
	double getOutput(double *pattern);
	double* getTransPos();
	double* getTransNeg();
	void calFLDA();
	void sampling();

	void transSamples(double **posSample, double **negSample, int nPos, int nNeg);

	double **centPos, **centNeg;
	double *sigPos, *sigNeg;
	double *weight;
	double J;
	int ndim, nCent;

	int nPos, nNeg;
	double **pos, **neg;
	double **posSample, **negSample;

	
	void getKernelizedSample(double **sample, int n, double **trans);
	void getFisherCriterion(double **pos, double **neg, int dim);
};

class CPONBDA{
public:
	CPONBDA();
	//Make Distributions by Samples
	CPONBDA(double *posSample, double *negSample, int nPos, int nNeg);

	//Make Distributions by Params
	CPONBDA(int nPos, int nNeg, int nPosSample, int nNegSample, double *posParam, double *negParam);
	~CPONBDA();

	double* getUncertaintyValue();
	double getUncertaintyMeasure();
	bool isUncertain(double val);
	double getLower();
	double getUpper();
	double getUncertainLower();
	double getUncertainUpper();
	void calUncertaintyValue();

	void printUncertainGraph(const char path[]);


	BD *pos, *neg;
	bool blp;
	
	double *uncertaintyValue;
	double uncertaintyMeasure;
	double getKa0(double x);
	double uncertainL, uncertainU;
	double L,U,theta;
};

class BD{
public:
	BD(double *sample, int n);
	BD(int n, int samplingSize, double *param);
	~BD();
	BD();

	double getCDFValue(double x);
	void getEstimator();
	int getSampleSize();
	double* getParameter();
	double getMean();

	double getKSPValue(double t);
	double dn, pValue;
	bool bSampling;
	int n, pn;
	double  mParam[4];
	double *sample;
	double *preSample;
	
	int SAMPLINGSIZE;

	// The sample data must be sorted.
	double* getDel(double *sample, int n, double *params);
	double getSdelsquared(double *sample, int n, double *params);
	double* getGradientSdel(double *sample, int n, double *params);
	double* getAlphakbyLS(double *sample, int n, double *params, double *dk);
	
	void getEstimatorbyCG();
	void getInitEstimator(double params[]);
	void getEstimatorbyPaper();

	// Sampling Data
	double* sampling();
	double* sampling(double per);
	
	// Kormogrov-Sminov Test
	double getKSStatistic();
};




/*
	Using simpson 3/8 Rule.
	It has limitation if a or b less than 1.
	PDF is not able to integrate.
	It's only for a and b greater than 1.

double getBetaPDF(double x, double params[]);
double getBetaCDF(double x, double params[]);
double getIntegrationSimpsonBeta(double(*f)(double x, double params[]), double x, double params[], int n);

////////////////////////Not Use//////////////////////////////


double getBetaPDF(double x, double params[])
{
	double betafunc = tgamma(params->a + params->b) / (tgamma(params->a) * tgamma(params->b));
	double bias = 1.E-4;
	double result = 0.f;
	//1 - x if x is very small, perhaps smaller than machine epsilon.

	if (x - params->L <= bias) x += bias;
	if (params->U - x <= bias) x -= bias;

	result = betafunc * pow(x - params->L, params->a - 1)*pow(params->U - x, params->b - 1) / pow(params->U - params->L, params->a + params->b - 1);

	return result;
}

double getBetaCDF(double x, double params[])
{
	int n;

	n = 999;
	return getIntegrationSimpsonBeta(&getBetaPDF, x, params, n);
}

double getIntegrationSimpsonBeta(double(*f)(double x, double params[]), double x, double params[], int n)
{
	int i;
	double x0, x1, x2, x3;
	double result = 0.f, h;

	h = (x - params->L) / (double)n;

	if (x == params->L) return 0.f;
	if (x == params->U) return 1.f;

	x3 = params->L;

	for (i = 0; i <= n - 3; i += 3){
		x0 = x3;
		x1 = x0 + h;
		x2 = x1 + h;
		x3 = x2 + h;
		result += 3.f * h / 8.f * (f(x0, params) + 3.f * f(x1, params) + 3.f * f(x2, params) + f(x3, params));
	}

	if (n % 3 != 0)
	{
		n = n % 3;
		x0 = x3;

		if (n == 1){
			x1 = x0 + h;
			result += h / 2.f * (f(x0, params) + f(x1, params));
		}
		else{ // n == 2
			x1 = x0 + h;
			x2 = x1 + h;
			result += h / 3.f * (f(x0, params) + 4.f * f(x1, params) + f(x2, params));
		}
	}

	return result;
}
*/

#endif