#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include "gamma.h"
using namespace std;
#define NR_N 25
#define AN_N 35
#pragma warning(disable:4996)


int main() {
	int i, j;
	vector<double> rsignals;
	vector<gpar> params;
	char *path = "bin/", *ctxt = ".txt", fdir[100];
	ifstream fin;
	ofstream fout("an_alsk.txt");
	double temp;

	setbuf(stdout, NULL);

	for (i = 34; i < AN_N; i++) {  //NR_N; i++) {
		sprintf(fdir, "%s%d%s", path, 200+i, ctxt);
		fin.open(fdir);
		if (!fin) continue;
		rsignals.clear();
		rsignals.resize(0);
		params.clear();
		params.resize(0);
		
		fin >> temp;
		while (!fin.eof()) {		
			fin >> temp;

			rsignals.push_back(temp);
		}
		fin.close();
		MLE_window(&params, 30, 10, rsignals);

		for (j = 0; j < params.size(); j++) {
			fout << params[j].a << ' ' << params[j].l << ' ' << params[j].s << ' ' << params[j].k << endl;
		}
	}
	return 0;
}