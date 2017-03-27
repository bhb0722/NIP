#pragma once
#ifndef GKFN_H
#define GKFN_H

#include <cstdio>
#include <cmath>

#define SGN(x) ((x) > 0 ? (1):(-1)) 
#define SQUARE(x) (if x)
#define pi 3.141596
#define id 30       /* max. dimension of input space */   
#define di 800      /* max. dimension of PFUs */
#define Nm 1300     /* max. number of time series data */      
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
	  int ���� ������
	  int get A()
	  void learn Gkfn()
	  void test()
  
*/

class GKFN {
public:
	~GKFN();
	GKFN() {}
	GKFN(char *filename, int E, int Tau, int PredictionStep, double TraningRate);
	void learn(int NumOfKernels, int NumOfEpochs, double errMargin = 1.f, double UBofSTD = 1.f);
	double getRsquared() {
		return rsq;
	}
	double getRMSE() {
		return rmse;
	}
	//void test();
	
private:

	void RECRUIT_FTN();
	void GENERAL_INVERSE();
	void PREDICTION();
	double OUTPUT(int i);


	void SIGMA_HN(int iter);
	void MEAN_HN(int iter);
	void OUTWEIGHT_HN(int iter);
	void INITIALIZE_SIGMA(int mode);
	
private:
	int  N, Nt, Nf, N1, N2, Np, Ns, P2, Td, Tk, Tp, St;
	double ic, si, K0;
	double *ts, **tsi, *tso; /* data */
	double *tse;
	double *s, *hn, *ek, *ck, **m, **Si; /* param  */
	double *o_sse;

	// evaluate
	double rsq;
	double rmse;
};


#endif


