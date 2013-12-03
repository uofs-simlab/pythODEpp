#ifndef BASE_METHOD_H
#define BASE_METHOD_H

#include <core/common.h>
#include <core/hash.h>
#include <core/paramvalue.h>
#include <core/vec.h>
#include <ivps/baseivp.h>

class BaseMethod {
protected:
	BaseIVP* _ivp;
	long* _acceptedSteps;
	FP* _dtOld;
	
	bool _accept;
	bool _sparse;
	bool _benchmark;
	
	FP _newtonFail;
	FP _newtonTol;

public:
	BaseMethod(Hash<ParamValue>& params, BaseIVP* ivp);
	virtual ~BaseMethod();

	void SetSolverVariables(long* acceptedSteps, FP* dtOld);
	void SetAccept(bool accept);
	bool Accept() const;

	virtual void PreStep(const FP tn, FP& dt, Vec<FP>& yn);
	virtual void PostStep(const FP tn, FP& dt, Vec<FP>& yn);
	virtual void Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew) = 0;
	virtual void UpdateTimestep();

	virtual FP CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol);
	
	virtual void GetStats(Hash<ParamValue>& params) const;

	virtual const char* GetName() const = 0;
	virtual long GetOrder() const = 0;
	virtual long GetAuxOrder() const;
};

#endif
