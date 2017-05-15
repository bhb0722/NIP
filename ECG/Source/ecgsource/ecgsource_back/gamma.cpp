#include "gamma.h"


void MLE_window(vector<gpar> *params, int win, int shift, vector<double> rsig) {
	int i, st, ed, n = rsig.size(), m;
	double val;
	vector<double> sample, rrint;
	gpar p;
	bool chk = true;
	st = 0;

	while (chk) {
		sample.clear();
		sample.resize(0);
		rrint.clear();
		rrint.resize(0);
		ed = st + win;

		for (i = 0; i < n; i++) {
			val = rsig[i];
			if (val >= st && val < ed) {
				sample.push_back(val);
				if (i == n - 1) chk = false;
			}
			else if (val > ed) break;
		}

		m = sample.size();

		for (i = 0; i < m - 1; i++) {
			val = sample[i + 1] - sample[i];

			if(val != 0.f)
				rrint.push_back(val); // RR interval
		}

		getParams(&p, rrint);
		p.s = getSkewness(rrint);
		p.k = getKurtosis(rrint);
		params->push_back(p);
		
 		st += shift;
	}
}

//MLE
void getParams(gpar *p, vector<double> rrint) {
	int n = rrint.size();
	int i;
	double s = getS(rrint), s3, sum = 0.f;
	double a, l;

	s3 = (s - 3.f) * (s - 3.f) + 24.f * s;

	a = (3.f - s + sqrt(s3)) / (12.f * s); // init a = k;

	double dela, inva, inva6, nu, de;

	for (i = 0; i < 250; i++) {
		inva = 1.f / a;
		inva6 = 1.f / (6.f * a + 1);
		nu = inva - 3.f*inva6 - s;
		de = 18.f * inva6*inva6 - inva*inva;
		dela = nu / de;
		a = a - dela;
	}

	for (i = 0; i < n; i++)
		sum += rrint[i];

	sum /= (double)n;

	l = a / sum;

	p->a = a;
	p->l = l;
}

double getS(vector<double>  rrint) {
	int n = rrint.size();
	int i;
	double sum = 0.f, logsum = 0.f;
	double s;

	for (i = 0; i < n; i++) {
		sum += rrint[i];
		logsum += log(rrint[i]);
	}

	sum /= (double)n;
	sum = log(sum);
	logsum /= (double)n;
	s = sum - logsum; // init

	return s;
}

double getDigamma(double k) {
	return k / 2.f - log(k);
}

double getSkewness(vector<double> sample) {
	int n = sample.size();
	int i;
	double m3 = 0.f, s = 0.f, m = 0.f, tmp;

	for (i = 0; i < n; i++)
		m += sample[i];

	m /= (double)n;

	for (i = 0; i < n; i++) {
		tmp = (sample[i] - m);
		m3 += tmp * tmp * tmp;
		s += tmp * tmp;
	}

	m3 /= (double)n;
	s /= (double)(n - 1);
	s = pow(s, 1.5f);

	return m3 / s;
}

double getKurtosis(vector<double> sample) {
	int n = sample.size();
	int i;
	double m4 = 0.f, m2 = 0.f, m = 0.f, tmp;

	for (i = 0; i < n; i++)
		m += sample[i];

	m /= (double)n;

	for (i = 0; i < n; i++) {
		tmp = (sample[i] - m) * (sample[i] - m);
		m4 += tmp * tmp;
		m2 += tmp;
	}
	
	m4 /= (double)n;
	m2 /= (double)n;
	m2 *= m2;

	return m4 / m2 - 3.f;
}