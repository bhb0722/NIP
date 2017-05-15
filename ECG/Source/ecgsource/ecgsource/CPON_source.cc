#include "CPON.h"

void CPON::train(int DIM, int npos, int nneg, double **pos, double **neg) {
	this->DIM = DIM;
	this->npos = npos;
	this->nneg = nneg;
	this->pos = pos;
	this->neg = neg;

	is_libsvm = false;

	clustering();
	kernelizing();
	lda();
	beta();
	getUncertainSample();

	is_train = true;
	//getMisclassfiedSample();
}

void CPON::getMisclassfiedSample() {
	int i, j, dim = 2; // # of cent = 2
	double val, *kval = (double*)malloc(sizeof(double) * 2);
	vector<int> tlab;

	//			posp[i] = innerP(dim, weight, posk[i]);

	for (i = 0; i < npos; i++) {
		getOutput(pos[i], val, kval);
		if (val < 0.5f) {
			tlab.push_back(i);
		}
	}

	mpos = tlab.size();

	poso = (double**)malloc(sizeof(double*) * mpos);

	for (i = 0; i < mpos; i++) {
		poso[i] = (double*)malloc(sizeof(double) * dim);

		for (j = 0; j < dim; j++)
			poso[i][j] = posk[tlab[i]][j];
	}

	tlab.clear();
	tlab.resize(0);

	for (i = 0; i < nneg; i++) {
		getOutput(neg[i], val, kval);
		val = 1.f - val;
		if (val < 0.5f) {
			tlab.push_back(i);
		}
	}

	mneg = tlab.size();

	nego = (double**)malloc(sizeof(double*) * mneg);

	for (i = 0; i < mneg; i++) {
		nego[i] = (double*)malloc(sizeof(double) * dim);

		for (j = 0; j < dim; j++)
			nego[i][j] = negk[tlab[i]][j];
	}

#ifdef PRINT_HIST
	ofstream fout("normal_uncertain.txt");

	for (i = 0; i < mpos; i++) {
		for (j = 0; j < dim; j++)
			fout << poso[i][j] << ' ';
		fout << endl;
	}

	fout.close();

	fout.open("abnormal_uncertain.txt");
	for (i = 0; i < mneg; i++) {
		for (j = 0; j < dim; j++)
			fout << nego[i][j] << ' ';
		fout << endl;
	}
	fout.close();
#endif
}

void CPON::getUncertainSample() {
	int i, j;
	double val;
	vector<int> tlab;

	mdim = cpos + cneg;
	maxo = (double*)calloc(mdim, sizeof(double));

	for (i = 0; i < npos; i++) {
		val = innerP(mdim, weight, posk[i]);
		if (bda->isUncertain(val)) {
			tlab.push_back(i);
		}
	}

	mpos = tlab.size();

	poso = (double**)malloc(sizeof(double*) * mpos);

	for (i = 0; i < mpos; i++) {
		poso[i] = (double*)malloc(sizeof(double) * mdim);

		for (j = 0; j < mdim; j++) {
			poso[i][j] = posk[tlab[i]][j];

			if (maxo[j] < poso[i][j])
				maxo[j] = poso[i][j];
		}
	}

	tlab.clear();
	tlab.resize(0);

	for (i = 0; i < nneg; i++) {
		val = innerP(mdim, weight, negk[i]);
		if (bda->isUncertain(val)) {
			tlab.push_back(i);
		}
	}

	mneg = tlab.size();

	nego = (double**)malloc(sizeof(double*) * mneg);

	for (i = 0; i < mneg; i++) {
		nego[i] = (double*)malloc(sizeof(double) * mdim);

		for (j = 0; j < mdim; j++) {
			nego[i][j] = negk[tlab[i]][j];

			if (maxo[j] < nego[i][j])
				maxo[j] = nego[i][j];
		}
	}

	/*for (i = 0; i < mpos; i++) {
		for (j = 0; j < mdim; j++)
			poso[i][j] /= maxo[j];
	}

	for (i = 0; i < mneg; i++) {
		for (j = 0; j < mdim; j++)
			nego[i][j] /= maxo[j];
	}*/

#ifdef PRINT_HIST
	ofstream fout("normal_uncertain.txt");

	for (i = 0; i < mpos; i++) {
		for (j = 0; j < mdim; j++)
			fout << poso[i][j] << ' ';
		fout << endl;
	}

	fout.close();

	fout.open("abnormal_uncertain.txt");
	for (i = 0; i < mneg; i++) {
		for (j = 0; j < mdim; j++)
			fout << nego[i][j] << ' ';
		fout << endl;
	}
	fout.close();
#endif

}

bool CPON::getOutput(double* pat, double &val, double* &kval) {
	int i;
	int kdim = cpos + cneg;
	double wval;
	double p, n;

	kval = (double*)malloc(sizeof(double) * kdim);

	for (i = 0; i < cpos; i++)
		kval[i] = getKerOutCov(pat, DIM, alpha, posm[i], pos_incov[i]);

	for (i = 0; i < cneg; i++)
		kval[i + cpos] = getKerOutCov(pat, DIM, alpha, negm[i], neg_incov[i]);

	wval = innerP(kdim, kval, weight);

	/*for (i = 0; i < kdim; i++)
		kval[i] /= maxo[i];*/

	p = posb->getCDFValue(wval);
	n = 1.f - negb->getCDFValue(wval);
	val = p / (p + n);

	return bda->isUncertain(wval);
}


double CPON::getKerOutCov(double *pat, int dim, double alpha, double *mean, double **incov) {
	int i, j;
	double sum = 0.f, flo = 1.e-10f;
	double *temp = (double*)malloc(sizeof(double) * dim);

	for (i = 0; i < dim; i++)
		temp[i] = 0.f;

	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++)
			temp[i] += (pat[j] - mean[j]) * incov[j][i];
	}

	for (i = 0; i<dim; i++)
		sum += temp[i] * (pat[i] - mean[i]);

	free(temp);

	sum *= -0.5f * alpha;
	sum = exp(sum);

	if (sum < flo)
		return flo;

	return sum;
}

double CPON::getKerOut(double *pat, int dim, double alpha, double *mean, double **var) {
	int i;
	double sum = 0.f, flo = 1.e-10f;

	for (i = 0; i < dim; i++)
		sum += (pat[i] - mean[i]) * (pat[i] - mean[i]) / var[i][i];

	sum *= -0.5f * alpha;
	sum = exp(sum);

	if (sum < flo)
		return flo;

	return sum;
}

int CPON::nearestCent(double **cent, int nCent, int dim, double *pat) {
	int i, j;
	double dist, mdist = -1, mi;

	for (i = 0; i < nCent; i++) {
		dist = 0.f;
		for (j = 0; j < dim; j++)
			dist += (cent[i][j] - pat[j]) * (cent[i][j] - pat[j]);

		if (mdist == -1 || mdist > dist) {
			mdist = dist;
			mi = i;
		}
	}

	return mi;
}

void CPON::getCentroid_LVQ(double **sample, double **result, double ***sigma, int n, int dim, int nCent, double learning_rate, int epoch) {
	int i, j, k = 0, mi;
	int *cnt = new int[nCent];// = (int*)malloc(sizeof(int) * nCent);

	for (i = 0; i < nCent; i++) {
		cnt[i] = 0;

		for (j = 0; j < dim; j++)
			result[i][j] = sample[i][j]; // w_i(0)
	}

	while (k++ < epoch) {
		for (i = 0; i < n; i++) {
			mi = nearestCent(result, nCent, dim, sample[i]);

			for (j = 0; j < dim; j++)
				result[mi][j] += (sample[i][j] - result[mi][j]) * learning_rate;
		}
	}

	for (i = 0; i < n; i++) {
		mi = nearestCent(result, nCent, dim, sample[i]);
		cnt[mi]++;

		for (j = 0; j < dim; j++) {
			for (k = 0; k < dim; k++) {
				sigma[mi][j][k] += (sample[i][j] - result[mi][j]) * (sample[i][k] - result[mi][k]);
			}
		}
	}

	for (i = 0; i < nCent; i++) {
		for (j = 0; j < dim; j++) {
			for (k = 0; k < dim; k++) {
				sigma[i][j][k] /= (double)cnt[i];
			}
		}
	}

	free(cnt);
}

void CPON::clustering() {
	int i, j;
	// void getCentroid(double **sample, int n, double **result, double **sigma, int* K, int ndim)

	posm = (double**)malloc(sizeof(double*) * cpos);
	negm = (double**)malloc(sizeof(double*) * cneg);
	posv = (double***)malloc(sizeof(double**) * cpos);
	negv = (double***)malloc(sizeof(double**) * cneg);
	pos_incov = (double***)malloc(sizeof(double**) * cpos);
	neg_incov = (double***)malloc(sizeof(double**) * cneg);

	for (i = 0; i < cpos; i++) {
		posm[i] = (double*)malloc(sizeof(double) * DIM);
		posv[i] = (double**)calloc(DIM, sizeof(double*));

		for (j = 0; j < DIM; j++)
			posv[i][j] = (double*)calloc(DIM, sizeof(double));
	}

	for (i = 0; i<cneg; i++) {
		negm[i] = (double*)malloc(sizeof(double) * DIM);
		negv[i] = (double**)calloc(DIM, sizeof(double*));

		for (j = 0; j < DIM; j++) {
			negv[i][j] = (double*)calloc(DIM, sizeof(double));
		}
	}

	//getCentroid(pos, npos, posm, posv, &cpos, DIM);

	if(cpos != 0)
		getCentroid_LVQ(pos, posm, posv, npos, DIM, cpos);
	//getCentroid(neg, nneg, negm, negv, &cneg, DIM);

	if(cneg != 0)
		getCentroid_LVQ(neg, negm, negv, nneg, DIM, cneg);

	//pos_incov = getIterativeMatrixInversion(posv,DIM);
	//neg_incov = getIterativeMatrixInversion(negv, DIM);
	for (i = 0; i < cpos; i++)
		pos_incov[i] = getIterativeMatrixInversion(posv[i], DIM);

	for (i = 0; i < cneg; i++)
		neg_incov[i] = getIterativeMatrixInversion(negv[i], DIM);
}

void CPON::kernelizing() { // diagonal terms of cov. matrix
	int i, j;

	posk = (double**)malloc(sizeof(double*) * npos);
	negk = (double**)malloc(sizeof(double*) * nneg);

	for (i = 0; i < npos; i++) {
		posk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

		for (j = 0; j < cpos; j++)
			posk[i][j] = getKerOutCov(pos[i], DIM, alpha, posm[j], pos_incov[j]);

		for (j = 0; j < cneg; j++)
			posk[i][j + cpos] = getKerOutCov(pos[i], DIM, alpha, negm[j], neg_incov[j]);
	}

	for (i = 0; i < nneg; i++) {
		negk[i] = (double*)malloc(sizeof(double) * (cpos + cneg));

		for (j = 0; j < cpos; j++)
		{
			if (j == 64)
				j = j;
			negk[i][j] = getKerOutCov(neg[i], DIM, alpha, posm[j], pos_incov[j]);

		}
		for (j = 0; j < cneg; j++)
			negk[i][j + cpos] = getKerOutCov(neg[i], DIM, alpha, negm[j], neg_incov[j]);
	}
}

void CPON::lda() {
	int i;
	int dim = cpos + cneg;
	FLDA flda;

	weight = (double*)malloc(sizeof(double) * dim);

	flda.weight = weight;
	flda.nPos = npos;
	flda.nNeg = nneg;
	flda.getFisherCriterion(posk, negk, dim);

	posp = (double*)malloc(sizeof(double) * npos);
	negp = (double*)malloc(sizeof(double) * nneg);

	for (i = 0; i < npos; i++)
		posp[i] = innerP(dim, weight, posk[i]);

	for (i = 0; i < nneg; i++)
		negp[i] = innerP(dim, weight, negk[i]);

#ifdef PRINT_HIST
	ofstream fout("normal_hist.txt");

	fout << npos << endl;
	for (i = 0; i < npos; i++)
		fout << posp[i] << endl;

	fout.close();

	fout.open("abnormal_hist.txt");
	fout << nneg << endl;
	for (i = 0; i < nneg; i++)
		fout << negp[i] << endl;
	fout.close();
#endif

	/*	for (int i = 0; i < dim; i++)
	cout << weight[i] << ' ';
	cout << endl;*/
}

void CPON::beta() {
	bda = new CPONBDA(posp, negp, npos, nneg);

	posb = bda->pos;
	negb = bda->neg;

#ifdef PRINT_PARAM
	ofstream fout("normal_param.txt");
	fout << posb->mParam[0] << endl;
	fout << posb->mParam[1] << endl;
	fout << posb->mParam[2] << endl;
	fout << posb->mParam[3] << endl;
	fout.close();
	fout.open("abnormal_param.txt");
	fout << negb->mParam[0] << endl;
	fout << negb->mParam[1] << endl;
	fout << negb->mParam[2] << endl;
	fout << negb->mParam[3] << endl;
	fout.close();
#endif
}