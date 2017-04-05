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
	instance : 클래스 내부 변수
	function : 함수
	     1) 생성자(Constructor) : 새로 객체를 만들때
		    GKFN *model = new GKFN();
		 2) 소멸자(Destructor) : 객체 내에서 사용된 인스턴스 변수들의 동적 할당을 해제 할 때 쓰임. => 메모리 누수 방지
		    free(model);
			delete model;
	     3) 연산자(Operator)
			ex) class MATRIX a, b
			    a+b
		 4) 임의 함수

	public : 외부 및 내부(상속) 클래스 모두 접근 가능
	protected : 외부애서는 불가. 내부(상속) 클래스 접근 가능
	private : 외부 및 내부(상속) 클래스 모두 접근 불가

	캡슐화(encapsulation) -> setter / getter
	private int a;
	
	public void set(int a) { this->a = a; }
	public int getA() { return this->a; }

	변수 네이밍(naming)
	  int mSampleSize;
	  GKFN *mGKFN;
	함수 네이밍
	  int 동사 명사형
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



