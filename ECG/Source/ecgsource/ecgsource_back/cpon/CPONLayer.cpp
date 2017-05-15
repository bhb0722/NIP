/*
Copyright 2015. Bae Hanbin(bhb0722@nate.com or bhb0722@skku.edu) all rights reserved.
Beta Functions for Deep Decision Network
*/

#include "betadist.h"

double* CPONLAYER::getCPONOutput(double* input){
	int i,j;
	double *output = (double*)malloc(sizeof(double) * nodeSize);
	double deSum = 0.f;

	for (i = 0; i < nodeSize; i++){
		CPONNODE *inode = mNode[i];

		if (!inode->checkUncertainty(input[i] / 1000.f))
			output[i] = inode->getCPONOutput(input[i] / 1000.f, inode->mInput);
		else{
			vector<int> *dep = inode->dep;
			int ndim = inode->ndim;
			double *pattern = (double*)malloc(sizeof(double)*ndim);
			

			for (j = 0; j < ndim; j++){
				pattern[j] = input[ dep->at(j) ];
			}

			output[i] = inode->getCPONOutput(inode->mFLDA[0]->getOutput(pattern), inode->mOutput);

			free(pattern);
		}

		deSum += output[i];
	}

	for (i = 0; i < nodeSize; i++)
		output[i] /= deSum;

	return output;
}

void CPONLAYER::trainInputClass(){
	int i, k;
	double p;
	double posParam[3347][4], negParam[3347][4];
	char c1[20] = "save/iposparams.txt",c2[20]="save/inegparams.txt",c3[20]="save/classprob.txt";
	string str = outDir + (string)c1;
	ifstream fin(str.c_str());
	str = outDir + (string)c2;
	ifstream fin2(str.c_str());
	LOG("trainInputClass", "");

	for (i = 0; i < nodeSize; i++){
		int l, nPos, nNeg;
		double dn;
		if (i % 1000 == 0) LOG("mInput",i);
		fin >> l >> nPos >> dn >> p;
		fin >> posParam[l][0] >> posParam[l][1] >> posParam[l][2] >> posParam[l][3];

		fin2 >> l >> nNeg >> dn >> p;
		fin2 >> negParam[l][0] >> negParam[l][1] >> negParam[l][2] >> negParam[l][3];
	}

	for (i = 0; i < nodeSize; i++){
		if (i== 43){
			mNode[i] = new CPONNODE(i, nodeSize, new CPONBDA(1000, 100, 1000, 100, posParam[i], negParam[i]));
		}
		else if (i == 0 || i == 42)
			mNode[i] = new CPONNODE(i, nodeSize, new CPONBDA(1000, 500, 1000, 500, posParam[i], negParam[i]));
		else
			mNode[i] = new CPONNODE(i, nodeSize, new CPONBDA(1000, 1000, 1000, 1000, posParam[i], negParam[i]));
	}

	fin.close();
	fin2.close();

	str = outDir + (string)c3;
	fin.open(str.c_str());

	for (i = 0; i < nodeSize; i++){
		fin >> p;
		mNode[i]->classProb = p;
	}

	fin.close();

}

void CPONLAYER::trainInputClass2(){
	int i, n, nPosSample, nNegSample;
	double L, U, theta;
	double posParam[4], negParam[4];
	string str = outDir + (string) "save/iparams.txt";
	ifstream fin(str.c_str());

	LOG("trainInputClass", "");

	for (i = 0; i < nodeSize; i++){

		if (i % 1000 == 0) LOG("mInput", i);


		fin >> n >> L >> U >> theta;
		fin >> nPosSample >> posParam[0] >> posParam[1] >> posParam[2] >> posParam[3];
		fin >> nNegSample >> negParam[0] >> negParam[1] >> negParam[2] >> negParam[3];
		mNode[i]->mInput = new CPONBDA(nPosSample, nNegSample, nPosSample, nNegSample, posParam, negParam);
		mNode[i]->mInput->calUncertaintyValue();
	}

	fin.close();

	str = outDir + (string)"save/classprob.txt";
	fin.open(str.c_str());

	for (i = 0; i < nodeSize; i++){
		double p;
		fin >> p;
		mNode[i]->classProb = p;
	}

	fin.close();
}

void CPONLAYER::trainOutputClass(){
	int i,n, nPosSample, nNegSample;
	double L, U, theta;
	double posParam[4], negParam[4];
	string str = outDir + (string)"save/oparams.txt";
	ifstream fin(str.c_str());

	LOG("trainOutputClass", "");

	for (i = 0; i < nodeSize; i++){

		if (i % 1000 == 0) LOG("mOutput", i);

		
		fin >> n >> L >> U >> theta;
		fin >> nPosSample >> posParam[0] >> posParam[1] >> posParam[2] >> posParam[3];
		fin >> nNegSample >> negParam[0] >> negParam[1] >> negParam[2] >> negParam[3];
		mNode[i]->mOutput = new CPONBDA(nPosSample, nNegSample, nPosSample, nNegSample, posParam, negParam);
		mNode[i]->mOutput->L = L;
		mNode[i]->mOutput->U = U;
		mNode[i]->mOutput->theta = theta;
		//mNode[i]->mOutput->calUncertaintyValue();
	}

	fin.close();
}

void CPONLAYER::trainFLDA(){
	int i,j,k,z;
	char t1[100] = "cdat/flda_", t2[20] = ".txt",
		t3[100] = "save/urelated.txt";

	stringstream sin;
	string str;
	ifstream fin;

	 LOG("trainFLDA", "");

	 str = outDir + (string)t3;

	 fin.open(str.c_str());

	 for (i = 0; i < nodeSize; i++)
	 {
		 CPONNODE *inode = mNode[i];
		 fin >> k;

		 for (j = 0; j < k; j++){
			 fin >> z;
			 inode->pushRelatedNode(z);
		 }


		 inode->mFLDA = (FLDA**)malloc(sizeof(FLDA*));
		 inode->mFLDA[0] = NULL;
	 }

	 fin.close();

	for(i=0;i<nodeSize;i++){

		if (i % 1000 == 0) LOG("mFLDA", i);

		sin.str("");
		sin.clear();
		sin << i;
		str = outDir + (string)t1 + sin.str() + (string)t2;
		fin.open(str.c_str());		

		CPONNODE *inode = mNode[i];
		int nCent, tdim, ndim;

		//LOG("trainFLDA", str.c_str());
		fin >> ndim >> nCent;
		//LOG(ndim, nCent);

		inode->ndim = ndim;
		inode->opti = 0;
		tdim = nCent * 2;



		double **centPos = (double**)malloc(sizeof(double*) * nCent);
		double **centNeg = (double**)malloc(sizeof(double*) * nCent);
		double *sigPos = (double*)malloc(sizeof(double) * nCent);
		double *sigNeg = (double*)malloc(sizeof(double) * nCent);
		double *weight = (double*)malloc(sizeof(double) * tdim);

		for (j = 0; j < nCent; j++)
			centPos[j] = (double*)malloc(sizeof(double) * ndim);
			
		for (j = 0; j < nCent; j++)
			centNeg[j] = (double*)malloc(sizeof(double) * ndim);

		for (j = 0; j < nCent; j++){
			for (k = 0; k < ndim; k++)
				fin >> centPos[j][k];
		}

		for (j = 0; j < nCent; j++){
			for (k = 0; k < ndim; k++)
				fin >> centNeg[j][k];
		}

		for (j = 0; j < nCent; j++)
			fin >> sigPos[j];

		for (j = 0; j < nCent; j++)
			fin >> sigNeg[j];

		for (j = 0; j < tdim; j++)
			fin >> weight[j];

		//LOG("ALLOCATION", "COMPLETE");
		inode->mFLDA[0] = new FLDA(ndim,nCent,centPos,centNeg,sigPos,sigNeg,weight);
		fin.close();
	}

}

CPONLAYER::CPONLAYER(int nodeSize, const char outDir[])
{
	int i,j;
//	ifstream fin("/home/ETRI/CPON/rawdata/save/classprob.txt");

	this->nodeSize = nodeSize;
	this->outDir = outDir;
	///vPos = new vector<double>[nodeSize];
	//vNeg = new vector<double>[nodeSize];
	mNode = new CPONNODE*[nodeSize];


	for (i = 0; i < nodeSize; i++){
	//	double p;
	//	fin >> p;
		mNode[i] = new CPONNODE(i);
	//	mNode[i]->classProb = p;
	}

//	fin.close();

//	bCheck = (bool*)calloc(sizeof(bool), nodeSize);

/*	fin.open("/home/ETRI/CPON/save/nNode.txt");

	for (i = 0; i < 91; i++){
		fin >> j;
		bCheck[j] = true;
	}

	fin.close();

	trainInputClass();
	trainOutputClass();
	trainFLDA();*/
}

CPONLAYER::CPONLAYER(int sampleSize, int nodeSize, int *label, int *nPos, double **sample){
	int i, j, l;

	this->nodeSize = nodeSize;
	this->sampleSize = sampleSize; 
	this->sample = sample;
	this->label = label;
	this->nPos = nPos;
	
	mNode = new CPONNODE*[nodeSize];
	posSamples = (double**)malloc(sizeof(double*)*nodeSize);
	negSamples = (double**)malloc(sizeof(double*)*nodeSize);

	for (i = 0; i < nodeSize; i++){
		posSamples[i] = (double*)malloc(sizeof(double) * nPos[i]);
		negSamples[i] = (double*)malloc(sizeof(double) * (sampleSize - nPos[i]));
		nPos[i] = 0;
	}
	
	for (i = 0; i < sampleSize; i++){
		l = label[i];
		for (j = 0; j < nodeSize; j++){
			if (j == l){
				posSamples[j][nPos[j]] = sample[i][j];
				nPos[j]++;
			}
			else
				negSamples[j][i - nPos[j]] = sample[i][j];
		}
	}
}




void CPONLAYER::saveRelated(){
	int i, j, l;
	ofstream fout("save/sRelated.txt");
	CPONNODE *inode;
	vector<int> *idep;
	fout << nodeSize << endl;

	for (i = 0; i < nodeSize; i++){
		inode = mNode[i];
		idep = inode->dep;
		l = idep->size();
		fout << l << ' ';

		for (j = 0; j < l; j++)
			fout << idep->at(j) << ' ';
		fout << endl;
	}

	fout.close();
}

void CPONLAYER::saveInput(){
	int i, nPos, nNeg;
	BD *ibd;
	double *sample;
	CPONBDA *ibda;
	CPONNODE *inode;

	ofstream fout3("save/sParams.txt");
	LOG("CPONLAYER", "saveInput()");

	for (i = 0; i < nodeSize; i++){
		inode = mNode[i];
		ibda = inode->mInput;
		ibd = ibda->pos;

		if (!ibd->bSampling){
			sample = ibd->sample;
			nPos = ibd->n;
		}
		else
		{
			sample = ibd->preSample;
			nPos = ibd->pn;
		}

		fout3 << i << endl;
		fout3 << ibd->mParam[0] << ' ' << ibd->mParam[1] << ' ' << ibd->mParam[2] << ' ' << ibd->mParam[3] << ' ' << endl;

		ibd = ibda->neg;
		sample = ibd->sample;
		nNeg = ibd->n;

		if (!ibd->bSampling){
			sample = ibd->sample;
			nNeg = ibd->n;
		}
		else
		{
			sample = ibd->preSample;
			nNeg = ibd->pn;
		}

		fout3 << ibd->mParam[0] << ' ' << ibd->mParam[1] << ' ' << ibd->mParam[2] << ' ' << ibd->mParam[3] << ' ' << endl;
		fout3 << nPos << ' ' << nNeg << endl;
	}
}

void CPONLAYER::saveData()
{
	int i, j, k, l, ndim, nCent;
	BD *ibd;
	CPONBDA *ibda;
	CPONNODE *inode;
	FLDA *iflda;

	ofstream fout("save.txt");
	fout << nodeSize;

	LOG("CPONLAYER", "saveData()");

	for (i = 0; i < nodeSize; i++){
		LOG(i, "th node");
		LOG("per", nodeSize);

		inode = mNode[i];

		fout << i << ' ' << inode->opti << endl;
		iflda = (inode->mFLDA)[inode->opti];


		// mInput
		ibda = inode->mInput;
		ibd = ibda->pos;
		fout << ibd->mParam[0] << ' ' << ibd->mParam[1] << ' ' << ibd->mParam[2] << ' ' << ibd->mParam[3] << ' ' << endl;
		fout << ibd->pValue << endl;
		ibd = ibda->neg;
		fout << ibd->mParam[0] << ' ' << ibd->mParam[1] << ' ' << ibd->mParam[2] << ' ' << ibd->mParam[3] << ' ' << endl;
		fout << ibd->pValue << endl;

		// mOutput
		ibda = inode->mOutput;
		ibd = ibda->pos;
		fout << ibd->mParam[0] << ' ' << ibd->mParam[1] << ' ' << ibd->mParam[2] << ' ' << ibd->mParam[3] << ' ' << endl;
		fout << ibd->pValue << endl;
		ibd = ibda->neg;
		fout << ibd->mParam[0] << ' ' << ibd->mParam[1] << ' ' << ibd->mParam[2] << ' ' << ibd->mParam[3] << ' ' << endl;
		fout << ibd->pValue << endl;

		// mFLDA

		l = inode->dep->size();

		fout << l << endl;

		for (j = 0; j < l; j++){
			fout << inode->dep->at(j) << ' ';
		}

		fout << endl;

		ndim = iflda->ndim;
		nCent = iflda->nCent;
		fout << ndim << ' ' << nCent << endl;

		for (j = 0; j < nCent; j++){
			for (k = 0; k < ndim; k++)
				fout << (iflda->centPos)[j][k] << ' ';
		}
		fout << endl;

		for (j = 0; j < nCent; j++){
			for (k = 0; k < ndim; k++)
				fout << (iflda->centNeg)[j][k] << ' ';
		}
		fout << endl;

		for (j = 0; j < nCent; j++)
			fout << (iflda->sigPos)[j] << ' ';
		fout << endl;

		for (j = 0; j < nCent; j++)
			fout << (iflda->sigNeg)[j] << ' ';
		fout << endl;

		for (j = 0; j < nCent * 2; j++)
			fout << (iflda->weight)[j] << ' ';
		fout << endl;
		fout << iflda->J << endl;

	}
}

void CPONLAYER::doTraining()
{
	int i, j, k, l, lab;
	int *arr = NULL;
	vector<int> *vt = new vector<int>();
	double *pattern;
	ofstream fout("save/trans.txt"), fout2("save/UM_BnA.txt"), fout3("save/ratio.txt");


	LOG("CPONLayer", "doTraining");

#if ISINPUT
	for (i = 0; i < nodeSize; i++){
		mNode[i] = new CPONNODE(i + 1, nodeSize, posSamples[i], negSamples[i], nPos[i], sampleSize - nPos[i]);
		LOG("mNode", i);
		LOG(mNode[i]->mInput->getUncertainLower(), mNode[i]->mInput->getUncertainUpper());
		LOG("Uncertainty Measure", mNode[i]->mInput->getUncertaintyMeasure());
	}
	saveInput();
#else
	/*DNN OUTPUT TEST*/
	int nPos, nNeg;
	double posParam[4], negParam[4];
	ifstream fin("save/sParams.txt");


	for (i = 0; i < nodeSize; i++){
		LOG("createNode", i);
		fin >> k;
		fin >> posParam[0] >> posParam[1] >> posParam[2] >> posParam[3];
		fin >> negParam[0] >> negParam[1] >> negParam[2] >> negParam[3];
		fin >> nPos >> nNeg;
		//mNode[i] = new CPONNODE(i + 1, nodeSize, new CPONBDA(nPos, nNeg, posParam, negParam));
	}


#endif

	for (i = 0; i < nodeSize; i++)
		mNode[i]->pushRelatedNode(i);

	for (i = 0; i < sampleSize; i++){
		lab = label[i];

		if (!mNode[lab]->checkUncertainty(sample[i][lab]))
			continue;

		LOG("check uncertainty of sample", i);

		for (j = 0; j < nodeSize; j++){
			if (mNode[j]->checkUncertainty(sample[i][j]))
				mNode[j]->pushRelatedNode(lab);
		}
	}

	for (i = 0; i < nodeSize; i++){
		mNode[i]->sortRelatedNode();
	}

	for (k = 0; k <= 10; k++)
	{
		vt = mNode[k]->getDep();
		l = vt->size();

		if (l > 0)
			arr = &vt->at(0);

		for (i = 0; i < sampleSize; i++)
		{
			lab = label[i];

			if (!mNode[k]->checkUncertainty(sample[i][k])) continue;

			pattern = (double*)malloc(sizeof(double)*l);

			for (j = 0; j < l; j++)
				pattern[j] = sample[i][arr[j]];

			if (k == lab)
				mNode[k]->pushUPos(pattern);
			else if (mNode[lab]->checkUncertainty(sample[i][lab]))
				mNode[k]->pushUNeg(pattern);
			else
				free(pattern);
		}

		mNode[k]->getOptFLDA(l);
		mNode[k]->setTransDist();
		CPONBDA *mInput = mNode[k]->mInput, *mOutput = mNode[k]->mOutput;

		cout << k << ' ' << mNode[k]->dep->size() << endl;
		cout << mNode[k]->uPosSample.size() << ' ' << mNode[k]->uNegSample.size() << endl;
		cout << mNode[k]->uPosSample.size() << ' ' << mNode[k]->uNegSample.size() << endl;
		cout << mOutput->pos->mParam[0] << ' ' << mOutput->pos->mParam[1] << ' ' << mOutput->pos->mParam[2] << ' ' << mOutput->pos->mParam[3] << ' ' << mOutput->pos->pValue << endl;
		cout << mOutput->neg->mParam[0] << ' ' << mOutput->neg->mParam[1] << ' ' << mOutput->neg->mParam[2] << ' ' << mOutput->neg->mParam[3] << ' ' << mOutput->neg->pValue << endl;
		cout << mInput->getUncertaintyMeasure() << ' ' << mOutput->getUncertaintyMeasure() * (mInput->getUncertainUpper() - mInput->getUncertainLower()) << endl;

		fout << k << ' ' << mNode[k]->dep->size() << endl;
		fout << mNode[k]->uPosSample.size() << ' ' << mNode[k]->uNegSample.size() << endl;
		fout << mOutput->pos->mParam[0] << ' ' << mOutput->pos->mParam[1] << ' ' << mOutput->pos->mParam[2] << ' ' << mOutput->pos->mParam[3] << ' ' << mOutput->pos->pValue << endl;
		fout << mOutput->neg->mParam[0] << ' ' << mOutput->neg->mParam[1] << ' ' << mOutput->neg->mParam[2] << ' ' << mOutput->neg->mParam[3] << ' ' << mOutput->neg->pValue << endl;
		fout2 << k << ' ' << mInput->getUncertaintyMeasure() << ' ' << mOutput->getUncertaintyMeasure() * (mInput->getUncertainUpper() - mInput->getUncertainLower()) << endl;
		double r = (double)mNode[k]->uPosSample.size() / (double)mNode[k]->mInput->pos->n;
		fout3 << i << ' ' << mNode[k]->mInput->pos->n << ' ' << mNode[k]->uPosSample.size() << ' ' << r << endl;

		mNode[k]->deleteUSamples();
	}
}
