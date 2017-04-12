/*
Copyright 2015. Bae Hanbin(bhb0722@nate.com or bhb0722@skku.edu) all rights reserved.
Beta Functions for Fisher Linear Discriminent Analysis
*/
#include "betadist.h"

FLDA::FLDA(double **posSample, double **negSample, int nPos, int nNeg, int ndim, int nCent)
{
	int i;
	int tdim = nCent * 2;

	this->posSample = posSample;
	this->negSample = negSample;
	this->nPos = nPos;
	this->nNeg = nNeg;
	this->ndim = ndim;
	this->nCent = nCent;

	centPos = (double**) malloc(sizeof(double*) * nCent);
	centNeg = (double**) malloc(sizeof(double*) * nCent);
	sigPos = (double*) malloc(sizeof(double) * nCent);
	sigNeg = (double*)malloc(sizeof(double) * nCent);
	weight = (double*)malloc(sizeof(double) * tdim);

	for (i = 0; i < nCent;i++)
		centPos[i] = (double*)malloc(sizeof(double) * ndim);

	for (i = 0; i < nCent; i++)
		centNeg[i] = (double*)malloc(sizeof(double) * ndim);

	pos = (double**)malloc(sizeof(double*) * nPos);
	neg = (double**)malloc(sizeof(double*) * nNeg);

	for (i = 0; i < nPos; i++)
		pos[i] = (double*)malloc(sizeof(double) * tdim);

	for (i = 0; i < nNeg; i++)
		neg[i] = (double*)malloc(sizeof(double) * tdim);
}

FLDA::FLDA(int ndim, int nCent, double **centPos, double **centNeg, double *sigPos, double *sigNeg, double *weight){
	this->ndim = ndim;
	this->nCent = nCent;
	this->centPos = centPos;
	this->centNeg = centNeg;
	this->sigPos = sigPos;
	this->sigNeg = sigNeg;
	this->weight = weight;
}

void FLDA::transSamples(double **posSample, double **negSample, int nPos, int nNeg)
{
	this->posSample = posSample;
	this->negSample = negSample;
	this->nPos = nPos;
	this->nNeg = nNeg;
	getKernelizedSample(posSample, nPos, pos);
	getKernelizedSample(negSample, nNeg, neg);
}

FLDA::~FLDA()
{
	/*int i;

	for (i = 0; i < nNeg; i++)
		free(neg[i]);

	for (i = 0; i < nPos; i++)
		free(pos[i]);

	for (i = 0; i < nCent; i++)
		free(centNeg[i]);

	for (i = 0; i < nCent; i++)
		free(centPos[i]);

	free(neg);
	free(pos);
	free(weight);
	free(sigNeg);
	free(sigPos);
	free(centNeg);
	free(centPos);*/
}

double* FLDA::getTransPos()
{
	int i, j, dim = nCent * 2;
	double *result = (double*) malloc(sizeof(double)*nPos);

	for (i = 0; i < nPos; i++)
	{
		result[i] = 0.f;
		for (j = 0; j < dim; j++){
			result[i] += weight[j] * pos[i][j];
			
			if (isnan(result[i])) {
				cout << i << ' ' << j;
				return NULL;
			}
		}
	}

	return result;
}

double* FLDA::getTransNeg()
{
	int i, j, dim = nCent * 2;
	double *result = (double*)malloc(sizeof(double)*nNeg);

	for (i = 0; i < nNeg; i++)
	{
		result[i] = 0.f;
		for (j = 0; j < dim; j++)
			result[i] += weight[j] * neg[i][j];

		if (isnan(result[i])) {
			cout << i << ' ' << j;
			return NULL;
		}
	}

	return result;
}

double FLDA::getOutput(double *pattern)
{
	int i, j, tCent = nCent * 2;
	double dist;
	double *trans = (double*)malloc(sizeof(double)*tCent);

	for (i = 0; i < nCent; i++){
		dist = 0.f;
		for (j = 0; j < ndim; j++){
			dist += (pattern[j] - centPos[i][j])*(pattern[j] - centPos[i][j]);
		}

		trans[i] = exp(-0.5f * dist / sigPos[i]);
	}

	for (i = 0; i < nCent; i++){
		dist = 0.f;
		for (j = 0; j < ndim; j++){
			dist += (pattern[j] - centNeg[i][j])*(pattern[j] - centNeg[i][j]);
		}

		trans[nCent + i] = exp(-0.5f * dist / sigNeg[i]);
	}

	dist = 0.f;

	for (i = 0; i < tCent; i++)
		dist += weight[i] * trans[i];

	free(trans);

	return dist;
}

double FLDA::getJ(){
	return J;
}

void FLDA::calFLDA()
{
	//getCentroid(posSample, nPos, centPos, sigPos, nCent, ndim);
	//getCentroid(negSample, nNeg, centNeg, sigNeg, nCent, ndim);
	getKernelizedSample(posSample, nPos, pos);
	getKernelizedSample(negSample, nNeg, neg);
	getFisherCriterion(pos, neg, nCent);
}

void FLDA::getFisherCriterion(double **pos, double **neg, int dim)
{
	int i, j, k;
//	int dim = nCent * 2;
	double *meanPos = (double*)malloc(sizeof(double)*dim);
	double *meanNeg = (double*)malloc(sizeof(double)*dim);
	double *meanDiff = (double*)malloc(sizeof(double)*dim); 
	double *diff = (double*)malloc(sizeof(double)*dim);
	double **S1 = (double**)malloc(sizeof(double*)*dim); 
	double **S2 = (double**)malloc(sizeof(double*)*dim); 
	double **invSw;
	double nu = 0.f, de = 0.f;

	for (i = 0; i < dim; i++){
		meanPos[i] = 0.f;
		meanNeg[i] = 0.f;
		S1[i] = (double*)malloc(sizeof(double)*dim);
		S2[i] = (double*)malloc(sizeof(double)*dim);
		for (j = 0; j < dim; j++)
		{
			S1[i][j] = 0.f;
			S2[i][j] = 0.f;
		}
	}

	for (i = 0; i < nPos; i++){
		for (j = 0; j < dim; j++)
			meanPos[j] += pos[i][j];
	}

	for (i = 0; i < nNeg; i++){
		for (j = 0; j < dim; j++)
			meanNeg[j] += neg[i][j];
	}

	for (i = 0; i < dim; i++){
		meanPos[i] /= (double)nPos;
		meanNeg[i] /= (double)nNeg;
		meanDiff[i] = meanNeg[i] - meanPos[i];
	}

	for (i = 0; i < nPos; i++){
		for (j = 0; j < dim; j++)
			diff[j] = pos[i][j] - meanPos[j];

		for (j = 0; j < dim; j++){
			for (k = 0; k < dim; k++)
				S1[j][k] += diff[j] * diff[k];
		}
	}

	for (i = 0; i < nNeg; i++){
		for (j = 0; j < dim; j++)
			diff[j] = neg[i][j] - meanNeg[j];

		for (j = 0; j < dim; j++){
			for (k = 0; k < dim; k++)
				S2[j][k] += diff[j] * diff[k];
		}
	}

	for (i = 0; i < dim; i++){
		weight[i] = 0.f;
		for (j = 0; j < dim; j++){
			S1[i][j] += S2[i][j];

			if(S1[i][j] < 1e-20f) S1[i][j] = 1e-20f;
		}
		
	}

	//if(dim < 3000)
		invSw = getIterativeMatrixInversion(S1, dim);
	//else
	//	invSw = getMatrixInverse(S1,dim);

	for (i = 0; i < dim; i++){
		for (j = 0; j < dim; j++)
			weight[i] += invSw[i][j] * meanDiff[j];
	}


	for (i = 0; i < dim; i++){
		meanPos[i] = 0.f;

		for (j = 0; j < dim; j++)
			meanPos[i] += meanDiff[j] * invSw[j][i];
	}

	J = 0.f;
	for (i = 0; i < dim; i++)
		J += meanPos[i] * meanDiff[i];

	//cout << nCent << ' ' << J << endl;

	free(meanPos);
	free(meanNeg);
	free(meanDiff);

	for (i = 0; i < dim; i++){
		free(S1[i]);
		free(S2[i]);
		free(invSw[i]);
	}

	free(S1);
	free(S2);
	free(invSw);
}

void FLDA::getKernelizedSample(double **sample, int n, double **trans)
{
	int i, j, k;
	double dist;

	for (i = 0; i < n; i++){
		for (k = 0; k < nCent; k++){
			dist = 0.f;

			for (j = 0; j < ndim; j++)
				dist += (sample[i][j] - centPos[k][j]) * (sample[i][j] - centPos[k][j]);

			trans[i][k] = exp(-0.5f * dist / sigPos[k]);
		}


		for (k = 0; k < nCent; k++){
			dist = 0.f;

			for (j = 0; j < ndim; j++)
				dist += (sample[i][j] - centNeg[k][j]) * (sample[i][j] - centNeg[k][j]);

			trans[i][k + nCent] = exp(-0.5f * dist / sigNeg[k]);
		}
	}

}

