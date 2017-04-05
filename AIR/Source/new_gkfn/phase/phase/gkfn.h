#pragma once
#ifndef GKFN_H
#define GKFN_H

#include <cstdio>
#include <cmath>

#define SGN(x) ((x) > 0 ? (1):(-1)) 
#define SQUARE(x) (if x)
#define pi 3.141596
#define id 200       /* max. dimension of input space */   
#define DIM 600      /* max. dimension of PFUs */
#define Nm 20000     /* max. number of time series data */      
#define sc 0.5      /* initial set-up coeff. of sigma */
#define er 0.9      /* error decrement rate */
#define Ip 10       /* number of iterations for parameter estimation */
#define cut 0.01

/*
	instance : Ŭ���� ���� ����
	function : �Լ�
	     1) ������(Constructor) : ���� ��ü�� ���鶧
		    GKFN *model = new GKFN();
		 2) �Ҹ���(Destructor) : ��ü ������ ���� �ν��Ͻ� �������� ���� �Ҵ��� ���� �� �� ����. => �޸� ���� ����
		    free(model);
			delete model;
	     3) ������(Operator)
			ex) class MATRIX a, b
			    a+b
		 4) ���� �Լ�

	public : �ܺ� �� ����(���) Ŭ���� ��� ���� ����
	protected : �ܺξּ��� �Ұ�. ����(���) Ŭ���� ���� ����
	private : �ܺ� �� ����(���) Ŭ���� ��� ���� �Ұ�

	ĸ��ȭ(encapsulation) -> setter / getter
	private int a;
	
	public void set(int a) { this->a = a; }
	public int getA() { return this->a; }

	���� ���̹�(naming)
	  int mSampleSize;
	  GKFN *mGKFN;
	�Լ� ���̹�
	  int ���� �����
	  int get A()
	  void learn Gkfn()
	  void test()
  
*/

class GKFN {
public:
	~GKFN();
	GKFN() {}
	GKFN(int N, double* ts, int E, int Tau, int PredictionStep, double TraningRate, int aq);
	GKFN(char *filename, int E, int Tau, int PredictionStep, double TraningRate);
	void learn(int NumOfKernels, int NumOfEpochs, double errMargin = 1.f, double UBofSTD = 1.f);
	double getTrainRsquared() {
		return trainRsq;
	}
	double getTestRsquared() {
		return testRsq;
	}
	double getTrainRMSE() {
		return trainRmse;
	}
	double getTestRMSE() {
		return testRmse;
	}
	double getYTrue() {
		return yTrue;
	}
	double getYEst() {
		return yEst;
	}
	//void test();
	
private:

	void RECRUIT_FTN();
	void GENERAL_INVERSE();
	void PREDICTION();


	void SIGMA_HN(int iter);
	void MEAN_HN(int iter);
	void OUTWEIGHT_HN(int iter);
	void INITIALIZE_SIGMA(int mode);
	
private:
	int  N, Nt, Nf, N1, N2, Np, Ns, P2, Td, Tk, Tp, St, AQ;
	double Rt;
	double ic, si, K0;
	double *ts, **tsi, *tso; /* data */
	double *tse;
	double *s, *hn, *ek, *ck, **m, **Si; /* param  */
	double *o_sse;

	// evaluate
	double trainRsq, testRsq;
	double trainRmse, testRmse;
	double yTrue, yEst;
};


#endif



