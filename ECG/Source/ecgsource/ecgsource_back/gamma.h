#pragma once
#include <vector>
#include <cmath>
using namespace std;

struct gpar {
	double a, l;
	double s, k;
};

void MLE_window(vector<gpar> *params, int win, int shift, vector<double> rsig);
double getS(vector<double> rrint);
void getParams(gpar *p, vector<double> rrint);
double getSkewness(vector<double> sample);
double getKurtosis(vector<double> sample);
double getDigamma(double k);
