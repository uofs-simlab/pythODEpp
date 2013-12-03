#ifndef ADDITIVE_EXPRK_METHOD_H
#define ADDITIVE_EXPRK_METHOD_H

#include <methods/rk.h>

class AdditiveExpRK : public BaseMethod {
protected:
	void ExpMtv(const BaseMat<FP>* M, const FP t, const Vec<FP>& v, Vec<FP>& expMtv);

	unsigned short _exponential;
	unsigned short _classical;

public:
	AdditiveExpRK(Hash<ParamValue>& params, BaseIVP* ivp);
	virtual ~AdditiveExpRK();
};

class DIRKCF1 : public AdditiveExpRK {
public:
	DIRKCF1(Hash<ParamValue>& params, BaseIVP* ivp);
	~DIRKCF1();

	virtual void Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew);

	virtual const char* GetName() const;
	virtual long GetOrder() const;
};

class DIRKCF2 : public AdditiveExpRK {
public:
	DIRKCF2(Hash<ParamValue>& params, BaseIVP* ivp);
	~DIRKCF2();

	virtual void Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew);

	virtual const char* GetName() const;
	virtual long GetOrder() const;
};

class DIRKCF3 : public AdditiveExpRK {
protected:
	bool NewtonSolve(BaseMat<FP>* jac, Vec<FP>& guess, FP t, FP dt, const Vec<FP>& constant);

public:
	DIRKCF3(Hash<ParamValue>& params, BaseIVP* ivp);
	~DIRKCF3();

	virtual void Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew);

	virtual const char* GetName() const;
	virtual long GetOrder() const;
};

class ERKCF2 : public AdditiveExpRK {
public:
	ERKCF2(Hash<ParamValue>& params, BaseIVP* ivp);
	~ERKCF2();

	virtual void Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew);

	virtual const char* GetName() const;
	virtual long GetOrder() const;
};

#endif

