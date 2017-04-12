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
}
// Nonlinear Conjugate Gradient Method

// Empirical Function del
// del[0] = E[del_i]
double* BD::getDel(double *sample, int n, double *params) {
	int i;
	double a = params[0], b = params[1], L = params[2], U = params[3];
	double yi, bi;
	double *pdf = (double*)malloc(sizeof(double) * n);
	double *del = (double*)malloc(sizeof(double) * n);

	bi = tgamma(a + b) / (tgamma(a) * tgamma(b));

	pdf[0] = 0.f;
	for (i = 1; i < n; i++) {
		yi = (sample[i] - L) / (U - L);
		pdf[i] = bi * pow(yi, a - 1) * pow(1 - yi, b - 1);
	}

	del[0] = 0.f;

	for (i = 1; i < n; i++) {
		yi = (sample[i] - sample[i - 1]) / (U - L);
		del[i] = yi / 2.f * (pdf[i] + pdf[i - 1]);

		del[0] += del[i];
	}

	del[0] /= n;
	free(pdf);
	return del;
}

double BD::getSdelsquared(double *sample, int n, double *params) {
	int i;
	double *del = getDel(sample, n, params);
	double sum = 0.f;

	for (i = 1; i < n; i++)
		sum += (del[i] - del[0]) * (del[i] - del[0]);

	sum /= n - 1;
	free(del);
	return sum;
}

double* BD::getGradientSdel(double *sample, int n, double *params) {
	int i;
	double a = params[0], b = params[1], L = params[2], U = params[3];
	double yi, yii;
	double bi = tgamma(a + b) / (tgamma(a) * tgamma(b));
	double *del = getDel(sample, n, params);
	double *A = (double*)malloc(sizeof(double) * n);
	double *B = (double*)malloc(sizeof(double) * n);
	double *grad = (double*)malloc(sizeof(double) * 2);

	A[0] = 0.f;
	B[0] = 0.f;

	for (i = 1; i < n; i++) {
		yi = (sample[i] - L) / (U - L);
		yii = (sample[i - 1] - L) / (U - L);
		if (yi < 0) yi = 0;
		else if (yi > 1.f) yi = 1.f;

		if (yii < 0) yi = 0;
		else if (yii > 1.f) yi = 1.f;

		A[i] = 0.f;
		B[i] = 0.f;


		if (i != n-1) {
			A[i] += pow(yi, a - 1) * pow(1 - yi, b - 1)*log(yi);
			B[i] += pow(yi, a - 1) * pow(1 - yi, b - 1)*log(1 - yi);
		}

		if (i != 1) {
			A[i] +=	pow(yii, a - 1) * pow(1 - yii, b - 1)*log(yii);
			B[i] += pow(yii, a - 1) * pow(1 - yii, b - 1)*log(1 - yii);
		}

		A[0] += A[i];
		B[0] += B[i];
	}

	A[0] /= n;
	B[0] /= n;

	grad[0] = 0.f;
	grad[1] = 0.f;

	for (i = 1; i < n; i++) {
		yi = (sample[i] - sample[i - 1]) / (U - L);
		if (yi < 0) yi = 0;
		else if (yi > 1.f) yi = 1.f;

		grad[0] += yi * bi * (del[i] - del[0]) * (A[i] - A[0]);
		grad[1] += yi * bi * (del[i] - del[0]) * (B[i] - B[0]);
	}

	grad[0] *= -1.f;// / (n - 1);
	grad[1] *= -1.f;// / (n - 1);

	free(A);
	free(B);
	free(del);

	return grad;
}

double* BD::getAlphakbyLS(double *sample, int n, double *params, double *dk) {
	int i;
	double p[3][4], S[3];
	double *ak = (double*)malloc(sizeof(double) * 2);

	ak[0] = 1.f;
	ak[1] = 1.f;

	for (i = 0; i < 4; i++) {
		p[0][i] = params[i];
		p[1][i] = params[i]; 
		p[2][i] = params[i];
	}

	while (1) {
		p[2][0] = p[0][0] + ak[0] * dk[0];
		p[2][1] = p[0][1] + ak[1] * dk[1];
		p[1][0] = p[0][0] + 0.5f * ak[0] * dk[0];
		p[1][1] = p[0][1] + 0.5f * ak[1] * dk[1];

		for (i = 0; i < 3; i++)
			S[i] = getSdelsquared(sample, n, p[i]);

		if (S[0] >= S[2] && S[2] >= S[1]) { // 1 > 3 > 2
			ak[0] *= 0.75f;
			ak[1] *= 0.75f;
			break;
		}else if (S[2] >= S[0] && S[0] >= S[1]) { // 3 > 1 > 2
			ak[0] *= 0.25f;
			ak[1] *= 0.25f;
			break;
		}
		else if (S[0] >= S[1] && S[1] >= S[2]) { // 1 > 2 > 3
			ak[0] *= 2.f;
			ak[1] *= 2.f;
		}
		else {
			ak[0] *= -1.f;
			ak[1] *= -1.f;
		}

		
	}

	return ak ;
}

void BD::getEstimatorbyPaper() {
	int i, k;
	double nu = 1.f;
	double p[4];
	double *gradk, *dk = (double*)malloc(sizeof(double) * 2);


	/*sample = sampling(0.90f);
	n = pn * 0.90f;
	getInitEstimator(mParam);
	*/
	/*mParam[0] = 2;
	mParam[1] = 2;
	mParam[2] = sample[0];
	mParam[3] = sample[n - 1];*/

	SAMPLINGSIZE = 500;
	sample = sampling();
	n = 500;
	sort(sample, sample + n);
	getInitEstimator(mParam);
	/*mParam[0] = 2;
	mParam[1] = 2;
	mParam[2] = sample[0];
	mParam[3] = sample[n - 1];*/


	for (k = 0; k < n; k++) {
		gradk = getGradientSdel(sample, n, mParam);

		if (mParam[0] + nu * gradk[0] > 1.f) mParam[0] += nu * gradk[0];
		if (mParam[1] + nu * gradk[1] > 1.f) mParam[1] += nu * gradk[1];
		dn = getKSStatistic();
		pValue = getKSPValue(sqrt((double)n) * dn); //pValue = getKSPValue(sqrt((double)pn) * dn);
	}
}

void BD::getEstimatorbyCG() {
	int i,k;
	double *gradk, *gradk2, *dk = (double*)malloc(sizeof(double) * 2);
	double *ak, *temp;
	double bk;


	sample = preSample;
	n = pn;
	getInitEstimator(mParam);

	// step 1
	gradk = getGradientSdel(sample, n, mParam);

	// step 2
	dk[0] = -1.f * gradk[0];
	dk[1] = -1.f * gradk[1];

	for (k = 0; k < n; k++) {
		// step 3
		ak = getAlphakbyLS(sample, n, mParam, dk); // Line Search 구현하다가 그만둠..

		// step 4
		mParam[0] += ak[0] * dk[0];
		mParam[1] += ak[1] * dk[1];
		
		dn = getKSStatistic();
		pValue = getKSPValue(sqrt((double)n) * dn); //pValue = getKSPValue(sqrt((double)pn) * dn);

		// step 5
		gradk2 = getGradientSdel(sample, n, mParam);
		bk = (gradk2[0] * gradk2[0] + gradk2[1] * gradk2[1]) / (gradk[0] * gradk[0] + gradk[1] * gradk[1]);
		temp = gradk;
		gradk = gradk2;
		free(temp);

		// step 6
		dk[0] = -1.f*gradk2[0] + bk * dk[0];
		dk[1] = -1.f*gradk2[1] + bk * dk[1];
	}
}

void BD::getEstimator()
{	
	int sampleSizes[8] = {0, 1000, 500, 100, 50, 30  };
	//int sampleSizes[8] = { 30, 50, 100, 500, 1000, 5000  }; // fast

	double samPer[5] = { 0.99f, 0.98f, 0.97f, 0.96f, 0.95f };
	int iterC = 0;
	int si;

	sampleSizes[0] = pn;

	for ( si = 0; si < 5; si++){
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

		//sample = sampling(samPer[si]);
		//n = SAMPLINGSIZE;
			int i, j, k, nfunk = 0;
			double p[5][4], y[5];
			p[0][0] = 0;
			getInitEstimator(p[0]);

			/*if (p[0][0] > 10.f || p[0][1] > 10.f) {
				p[0][0] /= 10.f;
				p[0][1] /= 10.f;
			}*/
			
			for (i = 1; i < 5; i++){
				for (j = 0; j < 4; j++)
					p[i][j] = p[0][j];
			}

			for (k = 0; k < 5; k++) {
				p[1][2] = (p[1][2] - p[1][3]) * 0.9f + p[1][3];

				p[2][3] = (p[2][3] - p[2][2]) * 0.9f + p[2][2];

				for (i = 0; i < 3; i++)
					y[i] = getDWLS(sample, n, p[i]);

				amoeba(sample, n, p, y, (double)1.0E-6, &nfunk);

				for (i = 1; i < 3; i++) {
					for (j = 0; j < 4; j++)
						p[i][j] = p[0][j];
				}
			}

			for (k = 0; k < 5; k++){
				if (p[1][1] < p[1][0])
					p[1][0] *= 1.5f;
				else
					p[1][1] *= 1.5f;

				if(p[2][1] < p[2][0])
					p[2][1] *= 0.5f;
				else
					p[2][0] *= 0.5f;
				
				//p[3][2] = (p[3][2] - p[3][3]) * 0.5f + p[3][3];

				//p[4][3] = (p[4][3] - p[4][2]) * 0.5f + p[4][2];

				for (i = 0; i < 3; i++)
					y[i] = getDWLS(sample, n, p[i]);

				amoeba(sample, n, p, y, (double)1.0E-6, &nfunk);

				for (i = 1; i < 3; i++){
					for (j = 0; j < 4; j++)
						p[i][j] = p[0][j];
				}
			}

			for (i = 0; i < 4; i++)
				mParam[i] = p[0][i];

			dn = getKSStatistic();
			pValue = getKSPValue(sqrt((double)n) * dn); //pValue = getKSPValue(sqrt((double)pn) * dn);
			if (bSampling) free(sample);
		} while (bSampling && pValue < 0.05f && iterC < 5);

		if (pValue >= 0.05f) break;
	}
	cout << pValue << endl;

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

double * BD::sampling(double per) {
	int i, nm = n - 1, off;
	double fn = (double)n;
	double d1;
	double *data;

	
	SAMPLINGSIZE = pn * per ;

	off = pn - SAMPLINGSIZE;

	data = (double*)malloc(sizeof(double) * SAMPLINGSIZE);
	
	for (i = 0; i < SAMPLINGSIZE; i++)
		data[i] = preSample[i + off];

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
	double a = mParam[0], b = mParam[1];
	double cdf = getBetaCDF(sample[0], mParam); // double cdf = getBetaCDF(preSample[0], mParam);
	double fn = (double)n;
	double D = 1 / fn - cdf, d1, d2;


	for (i = 0; i < n; i++) {  //for (i = 0; i < pn; i++){
		//if (i != nm && sample[i] == sample[i + 1])
		//	continue;

		//double cdf = getBetaCDF(preSample[i], mParam);
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
	double diff = (double)1e-20;
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
