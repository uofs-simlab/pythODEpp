#ifndef RKC_H
#define RKC_H

#include <methods/basemethod.h>

class RKC2 : public BaseMethod {
protected:
	Vec<FP> _F0;

	FP _eta;
	long _m;
	long _maxStages;
	Vec<FP> _spRadGuess;

	// Stats
	long _statMaxStages;
	long _statMinStages;
	long _statSteps;
	long _statStages;

	FP EstimateJacobianSpectralRadius(FP t, const Vec<FP>& y, Vec<FP>& guess, const Vec<FP>& yp, long split = 0);

public:
	RKC2(Hash<ParamValue>& params, BaseIVP* ivp);

	virtual void PreStep(const FP tn, FP& dt, Vec<FP>& yn);
	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew);
	virtual FP CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol);

	virtual void GetStats(Hash<ParamValue>& params) const;

	virtual const char* GetName() const;
	virtual long GetOrder() const;
	virtual long GetAuxOrder() const;

	static FP Cheb1(long n, FP x);
	static FP Cheb2(long n, FP x);
	static FP Cheb1p(long n, FP x);
	static FP Cheb1pp(long n, FP x);

	static FP Cheb1pRecursive(long n, FP Ujm1);
	static FP Cheb1ppRecursive(long n, FP Tj, FP Ujm1, FP x);
};

class RKC1 : public RKC2 {
public:
	RKC1(Hash<ParamValue>& params, BaseIVP* ivp);

	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew);

	virtual const char* GetName() const;
};

class PRKC : public RKC2 {
protected:
	Vec<FP> _Gm1;
	Vec<FP> _G0;
	Vec<FP> _Gmm1;
	Vec<FP> _Gm;
	Vec<FP> _K0;
	Vec<FP> _Km;
	Vec<FP> _Kf;

	FP _cmm1;
	Vec<FP> _spRadGuessG;
		
public:
	PRKC(Hash<ParamValue>& params, BaseIVP* ivp);
	
	virtual void PreStep(const FP tn, FP& dt, Vec<FP>& yn);
	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew);
	virtual FP CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol);
	
	virtual const char* GetName() const;
};

class IRKC : public RKC2 {
protected:
	FP _k1;
	Mat<FP>* _jac;

	void NewtonSolve(const Mat<FP>& LU, const Mat<FP>& P, const Vec<FP>& constant, FP t, FP dt, FP k1, Vec<FP>& k, Vec<FP>& Gj);

public:
	IRKC(Hash<ParamValue>& params, BaseIVP* ivp);

	virtual void PreStep(const FP tn, FP& dt, Vec<FP>& yn);
	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew);
	virtual FP CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol);
	
	virtual const char* GetName() const;
};

#endif

