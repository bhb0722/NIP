/*
	Copyright 2015. Bae Hanbin(bhb0722@nate.com or bhb0722@skku.edu) all rights reserved.
	Beta Functions for Beta Distribution
	*/

#include "betadist.h"

double BD::getCDFValue(double x){
	return getBetaCDF(x, mParam);
}

int BD::getSampleSize(){
	return SAMPLINGSIZE;
}

double* BD::getParameter(){
	return mParam;
}

BD::BD(int n, int samplingSize, double* param){
	mParam[0] = param[0];
	mParam[1] = param[1];
	mParam[2] = param[2];
	mParam[3] = param[3];
	this->n = n;
	this->preSample = NULL;
	this->SAMPLINGSIZE = samplingSize;
}

BD::BD(){
	n = 0;
	pn = 0;
	sample = NULL;
	preSample = NULL;
	bSampling = false;
}

BD::BD(double *sample, int n){
	//LOG("BD", "Constructor");
	this->preSample = sample;
	this->pn = n;
	sort(preSample, preSample + pn);
}

BD::~BD(){
	if (preSample)
		free(preSample);
}


void BD::getEstimator()
{	
	int sampleSizes[7] = {0, 5000, 1000, 800, 500, 300, 100};
	int iterC = 0;
	int si;

	sampleSizes[0] = pn;
	
	for ( si = 0; si < 7; si++){
		iterC = 0;
		SAMPLINGSIZE = sampleSizes[si];
		do{
			iterC++;
			if (pn == SAMPLINGSIZE){
				sample = preSample;
				n = pn;
				SAMPLINGSIZE = pn;
				bSampling = false;
			}
			else if(pn < SAMPLINGSIZE) break;
			else{
				sample = sampling();
				n = SAMPLINGSIZE;
				bSampling = true;
				sort(sample, sample + n);
			}
		
			int i, j, k, nfunk = 0;
			double p[5][4], y[5];
			p[0][0] = 0;
			getInitEstimator(p[0]);
			for (i = 1; i < 5; i++){
				for (j = 0; j < 4; j++)
					p[i][j] = p[0][j];
			}

			for (k = 0; k < 5; k++){
				p[1][0] += 1.f;

				p[2][1] += 1.f;

				p[3][0] += 0.5f;

				p[4][1] += 0.5f;

				for (i = 0; i < 5; i++)
					y[i] = getDWLS(sample, n, p[i]);

				amoeba(sample, n, p, y, (double)1.0E-6, &nfunk);

				for (i = 1; i < 5; i++){
					for (j = 0; j < 4; j++)
						p[i][j] = p[0][j];
				}
			}

			for (i = 0; i < 4; i++)
				mParam[i] = p[0][i];

			dn = getKSStatistic();
			pValue = getKSPValue(sqrt((double)n) * dn);
//			cout << pValue << endl;
			if (bSampling) free(sample);
		} while (bSampling && pValue < 0.05f && iterC < 20);

		if (pValue >= 0.05f) break;
	}
}



double* BD::sampling()
{
	//LOG("BD", "sampling()");
	int i;
	double *data;
	double rMax = (double)RAND_MAX;

	data = (double*)malloc(sizeof(double)*SAMPLINGSIZE);

	for (i = 0; i < SAMPLINGSIZE; i++){
		int j = (int)(rand() / rMax * (pn-1));
		data[i] = preSample[j];
	}

	return data;
}

void BD::getInitEstimator(double params[])
{
	//LOG("BD", "getInitEstimator()");
	int i;
	double sampleMean = 0.f, sampleVar = 0.f, L, U;
	double x, v;
	
	L = sample[0];
	U = sample[n - 1];

	for (i = 0; i < n; i++)
		sampleMean += sample[i];

	sampleMean = sampleMean / (double)n;

	for (i = 0; i < n; i++)
		sampleVar += (sample[i] - sampleMean)*(sample[i] - sampleMean);

	sampleVar = sampleVar / (double)(n - 1);

	//cout << sampleMean << ' ' << L << ' ' << U << endl;
	//cout << sampleVar << endl;

	x = (sampleMean - L) / (U - L);
	v = sampleVar / ((U - L) * (U - L));

	params[0] = x * (x * (1 - x) / v - 1) ;
	params[1] = (1 - x) * (x * (1 - x) / v - 1) ;

	if(params[0] < 0.f) params[0] *= -1.f;
	if(params[1] < 0.f) params[1] *= -1.f;

	params[2] = L;
	params[3] = U;
}


double BD::getMean(){
	double normMean = mParam[0] / (mParam[0] + mParam[1]);
	double region = mParam[3] - mParam[2];
	return mParam[3] + normMean * region;
}

double BD::getKSStatistic(){
	int i, nm = n - 1;
	//double a = params[0], b = params[1], L = params[2], U = params[3];
	double a = mParam[0], b = mParam[1];
	double cdf = getBetaCDF(sample[0], mParam);
	double fn = (double)n;
	double D = 1 / fn - cdf, d1, d2;


	for (i = 0; i < n; i++){
		if (i != nm && sample[i] == sample[i + 1])
			continue;

		double cdf = getBetaCDF(sample[i], mParam);

		d1 = fabs(cdf - (double)(i + 1) / fn);
		d2 = fabs(cdf - (double)i / fn);
		if (d1 > d2) d1 = d2;

		if (D < d1)
			D = d1;
	}

	return D;
}

double BD::getKSPValue(double t){
	int i = 1;
	double cur, sum = 0;
	double diff = (double)1e-100;
	double pi = 3.14159265359f;
	double coeff1 = (-1) * pi*pi / (8 * t*t);
	bool tf = true;

	do{
		cur = exp(coeff1 * (2 * i - 1) * (2 * i - 1));

		if (i > 1 && diff >= cur)
			tf = false;

		sum += cur;
		i++;
	} while (tf);

	if (t < FPMIN) t = FPMIN;

	sum *= sqrt(2 * pi) / t;

	return 1.0f - sum;
}
