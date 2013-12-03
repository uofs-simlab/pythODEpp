#ifndef ERK_METHOD_H
#define ERK_METHOD_H

#include <core/common.h>
#include <core/hash.h>
#include <core/paramvalue.h>
#include <core/vec.h>
#include <core/mat.h>
#include <methods/basemethod.h>

class RKMethod : public BaseMethod {
protected:
	Mat<FP> _a;
	Vec<FP> _b;
	Vec<FP> _baux;
	Vec<FP> _c;
	Vec<FP>* _k;
	
	long _m;
	
	void FillC();
	
public:
	RKMethod(Hash<ParamValue>& params, BaseIVP* ivp, long m);
	virtual ~RKMethod();

	const Mat<FP>& GetA() const;
	const Vec<FP>& GetB() const;
	const Vec<FP>& GetC() const;

	virtual FP CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol);
};

class ERK : public RKMethod {
public:
	ERK(Hash<ParamValue>& params, BaseIVP* ivp, long m);
	
	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew);
};

class DIRK : public RKMethod {
protected:
	void NewtonSolve(const FP tn, const FP dt, const long s, const Vec<FP>& yn, BaseMat<FP>* mat,
					 Vec<FP>& k, unsigned short split = 0);
	
public:
	DIRK(Hash<ParamValue>& params, BaseIVP* ivp, long m);
	
	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew);
};

class IMEX : public DIRK {
protected:
	Mat<FP> _a2;
	Vec<FP> _b2;
	Vec<FP> _baux2;
	Vec<FP> _c2;
	Vec<FP>* _k2;

	void FillC2();

public:
	IMEX(Hash<ParamValue>& params, BaseIVP* ivp, long m);
	virtual ~IMEX();

	const Mat<FP>& GetA2() const;
	const Vec<FP>& GetB2() const;
	const Vec<FP>& GetC2() const;
	
	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew);
	virtual FP CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol);
};

#endif

