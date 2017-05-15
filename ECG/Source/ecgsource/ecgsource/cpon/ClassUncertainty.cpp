/*
Copyright 2015. Bae Hanbin(bhb0722@nate.com or bhb0722@skku.edu) all rights reserved.
Beta Functions for CPONBDA
*/

#include "betadist.h"

CPONBDA::CPONBDA()
{
	pos = NULL;
	neg = NULL;
	uncertaintyValue = NULL;
}

CPONBDA::CPONBDA(int nPos, int nNeg, int nPosSample, int nNegSample, double *posParam, double *negParam){
	pos = new BD(nPos, nPosSample, posParam);
	neg = new BD(nNeg, nNegSample, negParam);
#if UV
	//calUncertaintyValue();
#endif
}

bool CPONBDA::isUncertain(double val) {
	return (uncertainL <= val && uncertainU >= val);
}

CPONBDA::CPONBDA(double *posSample, double *negSample, int nPos, int nNeg)
{
	uncertaintyValue = NULL;

	pos = new BD(posSample, nPos);
	neg = new BD(negSample, nNeg);

	pos->getEstimator();
	neg->getEstimator();

	calUncertaintyValue();

	//free(posSample);
	//free(negSample);
}

CPONBDA::~CPONBDA()
{
	//LOG("CPONBDA", "Destructor");
	delete pos;
	delete neg;

	delete[] uncertaintyValue;
}

void CPONBDA::printUncertainGraph(const char path[]){
	ofstream fout(path);

	double L = getLower(), U = getUpper();
	double h = (U - L) / (double)SIMPSONSTEP;

	fout << getUncertainLower() << ' ' << getUncertainUpper() << endl;

	for (int i = 0; i < SIMPSONSTEP; i++){
		fout << L << " " << uncertaintyValue[i] << endl;
		L += h;
	}

	fout.close();
}


double CPONBDA::getKa0(double x){
	double cdfPos = getBetaCDF(x, pos->getParameter());
	double cdfNeg = 1.f - getBetaCDF(x, neg->getParameter());
	double nPos = 1.f / sqrt((double)pos->n);
	double nNeg = 1.f / sqrt((double)neg->n);

	if (cdfPos > cdfNeg)
		return (cdfPos - cdfNeg) / (nPos + nNeg);
	else
		return (cdfNeg - cdfPos) / (nPos + nNeg);
}

double CPONBDA::getUncertaintyMeasure(){
	return uncertaintyMeasure;
}

double* CPONBDA::getUncertaintyValue(){
	return uncertaintyValue;
}

double CPONBDA::getLower(){
	return neg->mParam[2];
}

double CPONBDA::getUpper(){
	return pos->mParam[3];
}

double CPONBDA::getUncertainLower(){
	return uncertainL;
}

double CPONBDA::getUncertainUpper(){
	return uncertainU;
}

void CPONBDA::calUncertaintyValue()
{
	//LOG("CPONBDA", "getUncertaintyValue()");
	int i;
	int n = SIMPSONSTEP, nm = SIMPSONSTEP - 3;
	double x0, x1, x2, x3;
	double result = 0.f, h;
	double *posParam = pos->getParameter(), *negParam = neg->getParameter();
	double L = negParam[2], U = posParam[3], UL = posParam[2], UU = negParam[3];
	double max = -1.f;

	uncertaintyValue = (double*)malloc(sizeof(double)*(SIMPSONSTEP + 1));
	h = (U - L) / (double)n;
	x0 = x3 = L;
	uncertainL = L;
	uncertainU = U;

	for (i = 0; i <= n; i++){

		if (x0 >= UL && x0 <= UU)
			uncertaintyValue[i] = pos->getKSPValue(getKa0(x0)) / 2.f;
		else
			uncertaintyValue[i] = 0.f;

		if (uncertaintyValue[i] >= CONFLEVEL)
		{
			uncertainU = x0;

			if (uncertainL == L)
				uncertainL = x0;
		}

		if (uncertaintyValue[i] > max){
			max = uncertaintyValue[i];
			theta = x0;
		}

		x0 += h;
	}

	for (i = 0; i <= nm; i += 3){
		x0 = x3;
		x1 = x0 + h;
		x2 = x1 + h;
		x3 = x2 + h;
		result += (uncertaintyValue[i] + 3.f * uncertaintyValue[i + 1]
			+ 3.f * uncertaintyValue[i + 2] + uncertaintyValue[i + 3]);
	}

	result *= 3.f * h / 8.f;
	uncertaintyMeasure = result / (U - L);
}