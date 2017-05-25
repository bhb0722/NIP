#include "gkfn.h"
#include "utility.h"
#include <string.h>

#define PRTRST 1

// time deley

GKFN::~GKFN() {
	int i;

	for (i = 0; i < DIM; i++) {
		delete m[i];
		delete Si[i];
	}

	delete s; delete hn; delete  ek; delete  ck; delete  m; delete  Si; delete  o_sse;

	for (i = 0; i < Nm; i++)
		delete tsi[i];

	delete tso; delete  tsi; delete  tse;
}

GKFN::GKFN(int N, double* ts, int E, int Tau, int PredictionStep, double TraningRate, int aq, char *spot) {
	int i, j, k;
	double Rt;
	double X;
	AQ = aq;
	strcpy(SPOT, spot);
	tso = new double[Nm];
	tsi = new double*[Nm];
	tse = new double[Nm];

	for (i = 0; i < Nm; i++)
		tsi[i] = new double[id];

	s = new double[DIM];
	hn = new double[DIM];
	ek = new double[Nm];
	ck = new double[DIM];
	m = new double*[DIM];
	Si = new double*[DIM];
	o_sse = new double[DIM];

	for (i = 0; i < DIM; i++) {
		m[i] = new double[id];
		Si[i] = new double[DIM];
	}

	this->N = N;
	this->ts = ts;
	this->Td = Tau;
	this->Tk = E;
	this->Tp = PredictionStep;
	this->Rt = TraningRate;

	Ns = N - (E - 1) * Tau - 1;
	k = (E - 1) * Tau + 1;

	for (i = 1; i <= Ns; i++)
	{
		tso[i] = ts[k + Tp];

		for (j = 0; j < E; j++) {
			tsi[i][j + 1] = ts[k - j * Tau];
		}

		k++;
	}

	Rt = TraningRate;
	Nt = Ns*Rt;
	Nf = Ns - Nt;
	K0 = 1.f;
}

GKFN::GKFN(char *filename, int E, int Tau, int PredictionStep, double TraningRate) {
	int i,j,k;
	double X;
	FILE *fp;

	ts = new double[Nm];
	tso = new double[Nm];
	tsi = new double*[Nm];
	tse = new double[Nm];

	for (i = 0; i < Nm; i++)
		tsi[i] = new double[id];

	s = new double[DIM];
	hn = new double[DIM];
	ek = new double[Nm];
	ck = new double[DIM];
	m = new double*[DIM];
	Si = new double*[DIM];
	o_sse = new double[DIM];

	for (i = 0; i < DIM; i++) {
		m[i] = new double[id];
		Si[i] = new double[DIM];
	}


	/* generate training and test data */
	if ((fp = fopen(filename, "r")) == NULL)
		fprintf(stderr, "File Open Error : %s \n", filename);

	for (i = 1; fscanf(fp, "%lf", &ts[i]) == 1; ++i);

	fclose(fp);
	N = i - 1;



	this->Td = Tau;
	this->Tk = E;
	this->Rt = TraningRate;

	Ns = N - Td*(Tk - 1);
	for (i = 1; i <= Tk; ++i) {
		k = Td*(i - 1);
		for (j = 1; j <= Ns; ++j) {
			k = k + 1;
			tsi[j][i] = ts[k];
		}
	}

	this->Tp = PredictionStep;

	Ns = Ns - Tp;
	St = Td*(Tk - 1) + Tp;
	k = St;
	for (i = 1; i <= Ns; ++i) {
		k = k + 1;
		tso[i] = ts[k];
	}
	Rt = TraningRate;
	Nt = Ns*Rt;
	Nf = Ns - Nt;
	K0 = 1.f;
}

void GKFN::learn(int NumOfKernels, int NumOfEpochs, double errMargin, double UBofSTD) {
	this->ic = errMargin;
	this->si = UBofSTD;
	this->N1 = NumOfKernels;
	this->N2 = NumOfEpochs;

	RECRUIT_FTN(); /* initial set-up of PFN */
	GENERAL_INVERSE();
	PREDICTION();
}

void GKFN::PREDICTION() { 
	int i, j, n = Nt + Nf;
	double yavg, ytrue, yest;
	double nu, de;
	double yk, x;
	double *pka = new double[Np+1];
	FILE *res_train, *res_test;
	char rfname_train[50], rfname_test[50];

#if PRTRST
	sprintf(rfname_train, "result/%s_%d_%d_%d_%.2lf_%d_train.csv", SPOT, AQ, Tk, Td, Rt, N1);

	res_train = fopen(rfname_train, "w");

	sprintf(rfname_test, "result/%s_%d_%d_%d_%.2f_%d_test.csv", SPOT, AQ, Tk, Td, Rt, N1);

	res_test = fopen(rfname_test, "w");

	if (res_train == NULL) {
		printf("[ERROR] Can't open file %s\n\n", rfname_train);
		return;
	}
#endif
	trainRsq = 0.f;
	testRsq = 0.f;
	trainRmse = 0.f;
	testRmse = 0.f;
	yavg = 0.f;
	nu = 0.f;
	de = 0.f;

	for (i = 1; i <= n; i++) {
		for (j = 1; j <= Np; ++j) {
			x = s[j];
			pka[j] = kernx(Tk, x, tsi[i], m[j]);
		}

		yk = inner(Np, pka, ck);

		if (yk < 0)
			yk = 0.f;

		tse[i] = yk;
	}

	delete pka;

	for (i = 1; i <= Nt; i++)
		yavg += tso[i];

	yavg /= Nt;
#if PRTRST
	fprintf(res_train, "REAL,ESTIMATE\n");
#endif
	for (i = 1; i <= Nt; i++) {
		ytrue = tso[i];
		yest = tse[i];

		trainRmse += (ytrue - yest)*(ytrue - yest);
		de += (ytrue - yavg)*(ytrue - yavg);
#if PRTRST
		fprintf(res_train, "%lf,%lf\n", ytrue, yest);
#endif
	}
#if PRTRST
	fclose(res_train);
#endif
	nu = trainRmse;
	trainRsq = 1.f - nu / de;

	trainRmse /= Ns; // 여기 n-> Ns
	trainRmse = sqrt(trainRmse);

	yavg = 0.f;
	de = 0.f;
	nu = 0.f;

	for (i = Nt+1; i <= n; i++)
		yavg += tso[i];

	yavg /= Nf;
#if PRTRST
	fprintf(res_test, "REAL,ESTIMATE\n");
#endif
	for (i = Nt+1; i <= n; i++) {
		ytrue = tso[i];
		yest = tse[i];

		testRmse += (ytrue - yest)*(ytrue - yest);
		de += (ytrue - yavg)*(ytrue - yavg);
#if PRTRST
		fprintf(res_test, "%lf,%lf\n", ytrue, yest);
#endif
	}
#if PRTRST
	fclose(res_test);
#endif
	nu = testRmse;
	testRsq = 1.f - nu / de;

	testRmse /= Nf; // 여기 n -> Nf
	testRmse = sqrt(testRmse);

	
}

void GKFN::RECRUIT_FTN()
{
	int   a, b, g, i, j, k, l, n, p, q;
	double c, x, y, z, ec, dk, sk;
	double yd, yk, Ck;
	double *a1 = new double[id], *a2 = new double[id], *a3= new double[id];
	double *aka = new double[DIM], *akb = new double[DIM], *pka = new double[DIM], *pkb = new double[DIM], *Bka = new double[DIM], *Bkb = new double[DIM];
	double **Psa = new double*[DIM], **Psi = new double*[DIM], **dA = new double*[DIM];

	for (i = 0; i < DIM; i++) {
		Psa[i] = new double[DIM];
		Psi[i] = new double[DIM];
		dA[i] = new double[DIM];
	}


	x = 1000.; // min
	y = 0.;  // max
	for (j = 1; j <= Nt; ++j) {
		z = tso[j];
		if (z < x) {
			x = z;
			a = j;
		}
		if (z > y) {
			y = z;
			b = j;
		}
	}
	/*
	* assign input data
	*/
	for (j = 1; j <= Tk; ++j) {
		a1[j] = tsi[b][j];
		m[1][j] = a1[j];
		a2[j] = tsi[a][j];
		m[2][j] = a2[j];
		a3[j] = a1[j] - a2[j];
	}

	/*
	* Find maximum standard deviation
	*/
	x = sc*norm(Tk, a3);
	/*
	if (x > si)
	x = si;
	*/
	s[1] = x;
	s[2] = x;
	Psa[1][1] = 1.;
	Psi[1][1] = 1.;
	ck[1] = tso[b];
	i = 2;
	ec = ic;

	/* recruiting gkfs */
	while (i <= N1) {
		for (l = 1; l <= Nt; ++l) {
			/* input-sample */
			for (j = 1; j <= Tk; ++j) {
				a1[j] = tsi[l][j];
				yd = tso[l];
				m[i][j] = a1[j];
			}
			/* calculation of sigma */
			/*
			x = si;
			*/
			for (j = 1; j <= i - 1; ++j) {
				for (k = 1; k <= Tk; ++k) {
					a2[k] = m[j][k];
					a3[k] = a1[k] - a2[k];
				}
				y = sc*norm(Tk, a3);
				/*
				if (y < x) x = y;
				*/
			}
			/*
			s[i] = x;
			*/
			s[i] = y;
			/* update pka and pkb */
			for (j = 1; j <= i; ++j) {
				for (k = 1; k <= Tk; ++k) {
					a2[k] = m[j][k];
				}
				pkb[j] = kernx(Tk, x, a2, a1);
			}
			for (j = 1; j <= i; ++j) {
				for (k = 1; k <= Tk; ++k) {
					a2[k] = m[j][k];
				}
				x = s[j];
				pka[j] = kernx(Tk, x, a1, a2);
			}
			/* check error */
			yk = inner(i - 1, pka, ck);
			ek[l] = yd - yk;
			x = ek[l];
			if (abx(x) > ec) {
				/* recruitment of a PFU */
				/* update Psa */
				for (j = 1; j <= i - 1; ++j) {
					Psa[i][j] = pka[j];
					Psa[j][i] = pkb[j];
				}
				Psa[i][i] = pka[i];
				/* calculation of Psi */
				AX2(i - 1, i - 1, Psi, pka, aka);
				AX1(i - 1, i - 1, Psi, pkb, akb);
				sk = inner(i - 1, pkb, aka);
				dk = pka[i] - sk;
				outer(i - 1, akb, aka, dA);
				Ck = 1. / dk;
				scale(i - 1, Ck, dA, dA);
				Madd(i - 1, i - 1, Psi, dA, Psi);
				for (j = 1; j <= i - 1; ++j) {
					Bka[j] = -aka[j] / dk;
					Bkb[j] = -akb[j] / dk;
				}
				for (j = 1; j <= i - 1; ++j) {
					Psi[i][j] = Bka[j];
					Psi[j][i] = Bkb[j];
				}
				Psi[i][i] = Ck;
				/* calculation of ck */
				c = ek[l] / dk;
				for (j = 1; j <= i - 1; ++j)
					ck[j] -= c*akb[j];

				ck[i] = c;
				i++;
			}
			if (i>N1) break;
		}
		/* estimation of rms-error */
		z = 0.0;
		for (j = 1; j <= Nt; ++j)
			z += ek[j] * ek[j];

		z = sqrt(z / (double)Nt);
		ec *= er;
	}
	Np = i - 1;

	for (i = 0; i < DIM; i++) {
		delete Psa[i]; delete  Psi[i]; delete  dA[i];
	}

	delete a1; delete  a2; delete  a3; delete  aka; delete  akb; delete  pka; delete  pkb; delete  Bka; delete  Bkb; delete  Psa; delete  Psi; delete  dA;
}

void GKFN::GENERAL_INVERSE()
{
	register int q, p, l, j, k;
	int i;
	double x, y, z, sk, dk, Ck;
	double *aka = new double[DIM], *akb = new double[DIM], *ds = new double[DIM], **dA = new double*[DIM];

	for (i = 0; i < DIM; i++)
		dA[i] = new double[DIM];

	for (q = 1; q <= N2; ++q) {
		INITIALIZE_SIGMA(1);

		for (p = 1; p <= Ip; ++p) {
			for (l = 1; l <= Nt; ++l) {

				SIGMA_HN(l);
				/* update Si */
				AX1(Np, Np, Si, hn, aka);
				sk = inner(Np, hn, aka);
				dk = 1. + sk;
				Ck = 1. / dk;
				outer(Np, aka, aka, dA);
				scale(Np, Ck, dA, dA);
				Msub(Np, Np, Si, dA, Si);

				/* adjust sigma of pfn */
				AX1(Np, Np, Si, hn, ds);
				for (j = 1; j <= Np; ++j) {
					s[j] += ds[j] * ek[l];
					if (s[j] < 0)
						s[j] = -s[j];
				}
			} /* for-loop ; l */
		} /* for-loop ; p */

		  /* Mean update */
		INITIALIZE_SIGMA(2);
		for (p = 1; p <= Ip; ++p) {
			for (l = 1; l <= Nt; ++l) {

				MEAN_HN(l);

				/* update Si */
				AX1(Np*Tk, Np*Tk, Si, hn, aka);
				sk = inner(Np*Tk, hn, aka);
				dk = 1. + sk;
				Ck = 1. / dk;
				outer(Np*Tk, aka, aka, dA);
				scale(Np*Tk, Ck, dA, dA);
				Msub(Np*Tk, Np*Tk, Si, dA, Si);

				/* adjust sigma of pfn */
				AX1(Np*Tk, Np*Tk, Si, hn, ds);
				for (j = 1; j <= Np; ++j)
					for (k = 1; k <= Tk; ++k)
						m[j][k] += ds[Tk*(j - 1) + k] * ek[l];
			}
		}

		/* output weight update */
		/* initialization of Si */
		INITIALIZE_SIGMA(1);
		for (p = 1; p <= Ip; ++p) {
			for (l = 1; l <= Nt; ++l) {

				OUTWEIGHT_HN(l);

				/* update Psi */
				AX1(Np, Np, Si, hn, aka);
				AX2(Np, Np, Si, hn, akb);
				dk = 1. + inner(Np, hn, aka);
				Ck = 1. / dk;
				outer(Np, aka, akb, dA);
				scale(Np, Ck, dA, dA);
				Msub(Np, Np, Si, dA, Si);

				/* update ck */
				AX1(Np, Np, Si, hn, aka);
				for (j = 1; j <= Np; ++j)
					ck[j] += ek[l] * aka[j];
			}

			/* estimation of rms-error */
			z = 0.0;
			for (j = 1; j <= Nt; ++j)
				z += ek[j] * ek[j];

			z = sqrt(z / (double)Nt);

			/* added and changed by kmjung */
			o_sse[(q - 1)*Ip + p - 1] = z;
		}
	}

	for (i = 0; i < DIM; i++)
		delete dA[i];
	delete aka; delete  akb; delete  ds; delete  dA;
}

void GKFN::SIGMA_HN(int iter)
{
	register int j, k;
	double x, y, yd, yk;
	double *a1 = new double[id], *a2 = new double[id], *a3 = new double[DIM];

	

	/* input-sample */
	for (j = 1; j <= Tk; ++j)
		a1[j] = tsi[iter][j];
	yd = tso[iter];

	/* update pka */
	for (j = 1; j <= Np; ++j) {
		for (k = 1; k <= Tk; ++k)
			a2[k] = m[j][k];
		x = s[j];
		a3[j] = kernx(Tk, x, a1, a2);
	}

	/* check error */
	yk = inner(Np, a3, ck);
	ek[iter] = yd - yk;

	/* update hn */
	for (j = 1; j <= Np; ++j) {
		x = 1. / (2.*s[j] * s[j]);
		y = 0.0;
		for (k = 1; k <= Tk; ++k)
			y += (a1[k] - m[j][k])*(a1[k] - m[j][k]);
		hn[j] = ck[j] * y*exp(-x*y) / (s[j] * s[j] * s[j]);
	}
	delete a1; delete  a2; delete  a3;
}

void GKFN::MEAN_HN(int iter)
{
	register int j, k;
	double x, y, z, yd, yk;
	double *a1 = new double[id], *a2 = new double[id], *a3 = new double[DIM];

	


	/* input-sample */
	for (j = 1; j <= Tk; ++j)
		a1[j] = tsi[iter][j];
	yd = tso[iter];

	/* update pka */
	for (j = 1; j <= Np; ++j) {
		for (k = 1; k <= Tk; ++k)
			a2[k] = m[j][k];
		x = s[j];
		a3[j] = kernx(Tk, x, a1, a2);
	}

	/* check error */
	yk = inner(Np, a3, ck);
	ek[iter] = yd - yk;

	/* update hn */
	for (j = 1; j <= Np; ++j) {
		x = ck[j] / (s[j] * s[j]);

		y = 0.0;
		for (k = 1; k <= Tk; ++k)
			y += (a1[k] - m[j][k])*(a1[k] - m[j][k]);
		y /= (2.0*s[j] * s[j]);

		z = x*exp(-y);
		for (k = 1; k <= Tk; ++k)
			hn[Tk*(j - 1) + k] = z*(a1[k] - m[j][k]);
	}

	delete a1; delete  a2; delete  a3;
}

void GKFN::OUTWEIGHT_HN(int iter)
{
	register int j, k;
	double x, y, yd, yk;
	double *a1 = new double[id], *a2 = new double[id];

	/* input-sample */
	for (j = 1; j <= Tk; ++j)
		a1[j] = tsi[iter][j];
	yd = tso[iter];

	/* update pka */
	for (j = 1; j <= Np; ++j) {
		for (k = 1; k <= Tk; ++k)
			a2[k] = m[j][k];
		x = s[j];
		hn[j] = kernx(Tk, x, a1, a2);
	}

	/* check error */
	yk = inner(Np, hn, ck);
	ek[iter] = yd - yk;

	delete a1; delete  a2;
}

void GKFN::INITIALIZE_SIGMA(int mode) /* 1 = Sigma, Outweight ;; 2 = Mean */
{
	register int j, k;

	/* initialization of Si */
	if (mode == 1) {
		for (j = 1; j <= Np; ++j)
			for (k = 1; k <= Np; ++k) {
				if (j == k)
					Si[j][k] = K0;
				else
					Si[j][k] = 0.;
			}
	}
	else if (mode == 2) {
		for (j = 1; j <= Np*Tk; ++j)
			for (k = 1; k <= Np*Tk; ++k) {
				if (j == k)
					Si[j][k] = K0;
				else
					Si[j][k] = 0.;
			}
	}
}


//void GKFN::test()
//{
//	register int i, j, k;
//	double m1, m2, x, z, w, yd, yk, sum;
//	double a1[id], a2[id], pka[DIM];
//	double tmp1[1000];
//	/* evaluation on test-samples */
//
//	/* normalization term for training samples */
//	m1 = .0;
//	for (i = 1; i <= Nt; ++i) {
//		m1 = m1 + tso[i];
//	}
//	m1 = m1 / (double)Nt;
//	x = .0;
//
//	for (i = 1; i <= Nt; ++i) {
//		x = x + (tso[i] - m1)*(tso[i] - m1);
//	}
//	m1 = sqrt(x / ((double)Nt - 1.));
//	/* normalization term for test samples */
//	m2 = .0;
//	for (i = Nt + 1; i <= Nt + Nf; ++i) {
//		m2 = m2 + tso[i];
//	}
//	m2 = m2 / (double)Nf;
//
//	x = .0;
//	for (i = Nt + 1; i <= Nt + Nf; ++i) {
//		x = x + (tso[i] - m2)*(tso[i] - m2);
//	}
//	m2 = sqrt(x / ((double)Nf - 1.));
//
//	/* calculation of normalized rms error */
//	z = .0;
//
//	printf("yd yk ek\n");
//	for (i = 1; i <= Nt + Nf; ++i) {
//		/* input-sample */
//		for (j = 1; j <= Tk; ++j) {
//			a1[j] = tsi[i][j];
//			yd = tso[i];
//		}
//		/* update pka */
//		for (j = 1; j <= Np; ++j) {
//			for (k = 1; k <= Tk; ++k) {
//				a2[k] = m[j][k];
//			}
//			x = s[j];
//			pka[j] = kernx(Tk, x, a1, a2);
//		}
//		/* check error */
//		yk = inner(Np, pka, ck);
//
//		if (yk >= 0)
//			tmp1[i] = inner(Np, pka, ck);
//		else
//			tmp1[i] = 0;
//
//		ek[i] = yd - yk;
//		printf("%f %f %f\n", yd, yk, ek[i]);
//		z = z + ek[i] * ek[i];
//		if (i == Nt) {
//			w = sqrt(z / (double)(Nt - 1));
//			z = 0.0;
//		}
//	}
//
//	z = sqrt(z / (double)(Nf - 1));
//
//	/* write simulation results on output-files */
//
//
//
//
//	printf("mean vector ..\n");
//	for (i = 1; i <= Np; ++i) {
//		for (j = 1; j <= Tk; ++j) {
//			printf("%f ", m[i][j]);
//		}
//		printf("\n");
//	}
//
//	printf("sigma and ck \n");
//	for (i = 1; i <= Np; ++i) {
//		printf("%f %f\n", s[i], ck[i]);
//	}
//
//	sum = 0.0;
//	for (i = 1; i <= Nt + Nf; ++i)
//		sum += fabs(ek[i]);
//	sum /= (double)Nt + Nf;
//
//
//
//	//msd = metrics.mean_squared_error(true_num_pickups, predicted_num_pickups)
//	//m=metrics.r2_score(true_num_pickups, predicted_num_pickups)
//	printf("R^2 score\n");
//	int tmp0 = Nt + Nf; // Nf : # of test data
//
//	tmp1[0] = 0.0;
//
//	for (i = 1; i <= n; i++) {
//		nu += (tso[i] - tmp1[i]) * (tso[i] - tmp1[i]);
//		de += (tso[i] - tavg) * (tso[i] - tavg);
//	}
//
//	rsq = 1.f - nu / de;
//
//	printf("R^2 = %lf\n", rsq);
//
//
//	//return; 
//	for (i = 1; i <= tmp0; i++) {
//		var1 += (tso[i] - avg1)*(tso[i] - avg1);
//		var2 += (tmp1[i] - avg2)*(tmp1[i] - avg2);
//	}
//
//	var1 /= tmp0;
//	var2 /= tmp0;
//
//	for (i = 1; i <= tmp0; i++) {
//		for (j = 1; j <= tmp0; j++) {
//
//			cov += (tso[i] - avg1)*(tmp1[j] - avg2);
//		}
//
//	}
//
//	//printf("%f", cov);
//	cov = cov / (tmp0 - 1);
//
//	double S = ((cov)*(cov)) / (var1 * var2);
//	double S2 = cov / sqrt(var1*var2);
//	printf("S %f\n", S2);
//	printf("mena %f %f", avg1, avg2);
//	printf("var and cov %f  %f  %f\n", var1, var2, cov);
//
//
//
//
//
//
//
//	printf("Mean of absolute error = %10.4f\n", sum);
//
//	printf("Simulation Results of a Network with GKFs\n");
//	printf("adaptation of sigmas, mean vectors and output weights\n");
//	printf("\n");
//	printf("time delay = %d\n", Td);
//	printf("\n");
//	printf("input dimension = %d\n", Tk);
//	printf("\n");
//	printf("forward prediction step = %d\n", Tp);
//	printf("\n");
//	printf("number of total time series data = %d\n", N);
//	printf("\n");
//	printf("number of training data = %d\n", Nt);
//	printf("\n");
//	printf("number of test data = %d\n", Nf);
//	printf("\n");
//	printf("initial error margin = %f\n", ic);
//	printf("\n");
//	printf("upper-bound of standard deviation = %f\n", si);
//	printf("\n");
//
//	/*
//	printf ( "initial value of diagonal element in Si = %f\n", K0);
//	printf ( "\n");*/
//	printf("number of epochs for recruiting GKFs = %d\n", N1);
//	printf("\n");
//	printf("number of epochs for parameter estimation = %d\n", N2);
//	printf("\n");
//	printf("number of potential functions = %d\n", Np);
//	printf("\n");
//	/*
//	printf ( "quality of inversion: %f %f\n", q0, q1);
//	printf ( "\n");*/
//	printf("rms error for trainig-samples = %f\n", w);
//	printf("\n");
//	printf("standard deviation of trainig-samples = %f\n", m1);
//	printf("\n");
//	x = (w / m1)*100.;
//	printf("normlized rms error for trainig-samples = %f percent\n", x);
//	printf("\n");
//	printf("rms error for test-samples = %f\n", z);
//	printf("\n");
//	printf("standard deviation of test-samples = %f\n", m2);
//	printf("\n");
//	x = (z / m2)*100.;
//	printf("normlized rms error for test-samples = %f percent\n", x);
//	printf("\n");
//}
