/*
Copyright 2015. Bae Hanbin(bhb0722@nate.com or bhb0722@skku.edu) all rights reserved.
Beta Functions for CPONNODE
*/

#include "betadist.h"

CPONNODE::CPONNODE(int iLabel){
	this->iLabel = iLabel;
	dep = new vector<int>();
	dep->clear();
	bDep = (bool*)calloc(3347, sizeof(bool));
	opti = -1;
}
CPONNODE::CPONNODE(int iLabel, int nodeSize, CPONBDA *input)
{
	this->iLabel = iLabel;
	mInput = input;
	dep = new vector<int>();
	dep->clear();
	bDep = (bool*)calloc(nodeSize, sizeof(bool));
}

CPONNODE::CPONNODE(int iLabel, int nodeSize, double *posSample, double *negSample, int nPos, int nNeg){
	LOG("CPONNODE", "Constructor");

	this->iLabel = iLabel;
	mInput = new CPONBDA(posSample, negSample, nPos, nNeg);
	dep = new vector<int>();
	dep->clear();
	bDep = (bool*)calloc(nodeSize, sizeof(bool));
}
double CPONNODE::getDNNOutputExp2(double x){
	double y, yy;
	double pattern[2];
	double ret;
	pattern[0] = mInput->pos->getCDFValue(x);
	pattern[1] = 1 - mInput->neg->getCDFValue(x);
	y = mFLDA[0]->getOutput(pattern);
	yy = (y - mOutput->theta) / (mOutput->getUpper() - mOutput->getLower());

	ret = exp(10*yy);

	return ret;
}

double CPONNODE::getDNNOutputCPON(double x){
	double y;
	double pattern[2];
	pattern[0] = mInput->pos->getCDFValue(x);
	pattern[1] = 1 - mInput->neg->getCDFValue(x);
	y = mFLDA[0]->getOutput(pattern);

	double fPos = mOutput->pos->getCDFValue(y);
	double fNeg = mOutput->neg->getCDFValue(y);
	double denom;

	if (fPos == 0) return 0;
	denom = (fPos  + (1 - fNeg) );

	if (denom < 1e-20) denom = 1e-20;

	return fPos  / denom;
}

double CPONNODE::getDNNOutputWeightedCPON(double x){
	double y;
	double pattern[2];
	pattern[0] = mInput->pos->getCDFValue(x);
	pattern[1] = 1 - mInput->neg->getCDFValue(x);
	y = mFLDA[0]->getOutput(pattern);
		
	double fPos = mOutput->pos->getCDFValue(y);
	double fNeg = mOutput->neg->getCDFValue(y);
	double denom;

	if (fPos == 0) return 0;
	denom = (fPos * classProb + (1 - fNeg) *(1 - classProb));

	if (denom < 1e-20) denom = 1e-20;

	return fPos * classProb  / denom;
}
double CPONNODE::getOriginalCPONOutput(double x, CPONBDA *ibda){
	double fPos = ibda->pos->getCDFValue(x);
	double fNeg = ibda->neg->getCDFValue(x);
	double denom;

	if (fPos < 1e-20) return 0;
	if (fNeg < 1e-20) fNeg = 0;


	denom = fPos + (1.f - fNeg) ;

	if (denom < 1e-20) denom = 1e-20;

	return fPos * classProb / denom;
}


double CPONNODE::getCPONOutput(double x, CPONBDA *ibda){
	double fPos = ibda->pos->getCDFValue(x);
	double fNeg = ibda->neg->getCDFValue(x);
	double denom;
	
	if (fPos < 1e-20) return 0;
	if (fNeg < 1e-20) fNeg = 0;
	

	denom = (fPos*classProb + (1.f - fNeg) * (1.f - classProb));
	
	if (denom < 1e-20) denom = 1e-20;

	return fPos * classProb / denom;
}

double CPONNODE::getLogCPONOutput(double x, CPONBDA *ibda){
	double fPos = ibda->pos->getCDFValue(x);
	double fNeg = ibda->neg->getCDFValue(x);
	double denom;

	if (fPos == 0) return 0;
	denom = log ( fPos*classProb + (1 - fNeg) * (1 - classProb) + 1.0e-20f );

	return log(fPos) + log(classProb) - denom;
}

void CPONNODE::setClassProb(double prob){
	this->classProb = prob;
}

void CPONNODE::deleteUSamples(){
	int i,l;
	double *sample;

	l = uPosSample.size();

	for (i = 0; i < l; i++){
		sample = uPosSample[i];
		free(sample);
	}

	l = uNegSample.size();

	for (i = 0; i < l; i++){
		sample = uNegSample[i];
		free(sample);
	}

	free(uTransPos);
	free(uTransNeg);
}

void CPONNODE::pushUPos(double *pattern){
	uPosSample.push_back(pattern);
}

void CPONNODE::pushUNeg(double *pattern){
	uNegSample.push_back(pattern);
}

vector<int>* CPONNODE::getDep(){
	return dep;
}

bool CPONNODE::checkUncertainty(double x)
{
	double L = mInput->getUncertainLower(), U = mInput->getUncertainUpper();

	if (L <= x && x <= U) return true;

	return false;
}

void CPONNODE::pushRelatedNode(int p)
{
	if (!bDep[p]){
		dep->push_back(p);
		bDep[p] = true;
	}
}

void CPONNODE::sortRelatedNode(){
	sort(dep->begin(), dep->end());
}

CPONNODE::~CPONNODE(){
}


void CPONNODE::getOptFLDA(int ndim){
	int i, ini = -1, l = CENTMAX, cmin = CENTMIN, j;
	int nPos, nNeg;
	double pMax = 0.f,umMin = 1.1f;
	double *tPos, *tNeg;

	this->ndim = ndim;
	nPos = uPosSample.size();
	nNeg = uNegSample.size();
	
	LOG("getOptFLDA", "optCent");
	LOG("nPos", nPos);
	LOG("nNeg", nNeg);
	LOG("ndim", ndim);

	if (l >= nPos)
		l = nPos - 1;

	if (l >= nNeg)
		l = nNeg - 1;

	if (l <= cmin)
		cmin = l;

	l -= cmin - 1;

	mFLDA = (FLDA**)malloc(sizeof(FLDA*)*l);

	this->ndim = ndim;


	for (int i = 0; i < l; i ++){
		mFLDA[i] = new FLDA(&uPosSample[0], &uNegSample[0], nPos, nNeg, ndim, i + cmin);
		mFLDA[i]->calFLDA();
		
		tPos = mFLDA[i]->getTransPos();
		tNeg = mFLDA[i]->getTransNeg();

		if (tPos == NULL || tNeg == NULL){
			i--;
			break;
			cout << "NULL" << endl;
		}
		int ct = 0;
		double min, uv = 1.1f;
		CPONBDA *ibda = NULL;

		do{
			CPONBDA *bda = new CPONBDA(tPos, tNeg, uPosSample.size(), uNegSample.size());

			min = bda->pos->pValue;
			if (min > bda->neg->pValue) min = bda->neg->pValue;

			double uu = bda->uncertaintyMeasure;

			if (min >= 0.05f &&  uu < uv){
				uv = uu;

				//if (!ibda) delete ibda;
				ibda = bda;
				cout << uu << ' ';
				break;
			}
			//else
				//delete bda;
			
			//cout << ct << ' ';
		} while (++ct < 3);

		if (!ibda) continue;

		cout << endl;
		cout << "# of cent : "  << (int)(i + 1) << endl;
		cout << "pPos : " << ibda->pos->pValue << ' ';
		cout << "pNeg : " << ibda->neg->pValue << ' ';
		cout << "UM : " << ibda->uncertaintyMeasure << endl;


		if (min >= 0.05f){
			if (umMin > uv){
				if (ini != -1){
					delete mFLDA[ini];
					free(uTransPos);
					free(uTransNeg);
				}

				umMin = uv;
				uTransPos = tPos;
				uTransNeg = tNeg;
				mOutput = ibda;
				ini = i;

				if (umMin == 0.f){
					cout << "UM 0 break" << endl;
					break;
				}

				//break;
			}
			else {
				delete mFLDA[i];
				free(tPos);
				free(tNeg);
				//break;
			}
		}
		else{
			delete mFLDA[i];
			free(tPos);
			free(tNeg);
		}
	}

	opti = ini;
	LOG("Optimal Cent", ini + cmin);
}

void CPONNODE::setTransDist()
{
	//LOG("CPONNODE", "setTransDist() : create mOutput");
	uTransPos = mFLDA[opti]->getTransPos();
	uTransNeg = mFLDA[opti]->getTransNeg();
	mOutput = new CPONBDA(uTransPos, uTransNeg, uPosSample.size(), uNegSample.size());
	//mOutput->calUncertaintyValue();
}

void CPONNODE::getTransDist(){
	mFLDA[opti]->transSamples(&uPosSample[0], &uNegSample[0], uPosSample.size(), uNegSample.size());
	uTransPos = mFLDA[opti]->getTransPos();
	uTransNeg = mFLDA[opti]->getTransNeg();
	mOutput = new CPONBDA(uTransPos, uTransNeg, uPosSample.size(), uNegSample.size());

}

double CPONNODE::getOutput(double *pattern)
{
	if (bCertain) return pattern[0];
	else return mFLDA[opti]->getOutput(pattern);
}
