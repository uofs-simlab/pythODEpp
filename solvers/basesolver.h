#ifndef BASE_SOLVER_H
#define BASE_SOLVER_H

#include <core/common.h>
#include <core/exception.h>
#include <core/hash.h>
#include <core/paramvalue.h>
#include <core/timer.h>
#include <ivps/baseivp.h>
#include <methods/basemethod.h>

class BaseSolver {
protected:
	FP _dt;
	FP _tn;
	FP _tf;
	FP _stretch;
	FP _dtOld;

	bool _printTime;
	bool _complete;
	FP _lastWriteTime;
	FP _minWriteTime;
	
	Vec<FP> _yn;
	Vec<FP> _ynew;
	std::string _outPath;
	
	BaseIVP* _ivp;
	BaseMethod* _method;

	long _maxSteps;
	long _steps;
	long _acceptedSteps;
	long _rejectedSteps;
	Timer _timer;
	
	void CheckMaxSteps();
	void WriteFile(long f);

	virtual void UpdateTimestep();

	void StartTiming();
	void EndTiming();
	
public:
	BaseSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp);
	virtual ~BaseSolver();
	
	virtual void DumpRunInfo(Hash<ParamValue>& params);
	virtual void RunSimulation() = 0;
	virtual const char* GetName() const = 0;
};

class StepControlSolver : public BaseSolver {
public:
	enum SCType {
		STANDARD = 0,
		PREDICTIVE = 1,
		MIXED = 2
	};

protected:
	bool _rejectStep;

	FP _rTol;
	FP _aTol;
	FP _safety;
	FP _minChange;
	FP _maxChange;
	FP _maxRejectedChange;

	FP _eps;
	FP _epsLast;
	SCType _scType;
	bool _restrictReject;

	FP StandardStepControl();
	FP PredictiveStepControl();

	FP CalculateFactor();
	
	bool CheckMethodReject();
	void AcceptStep();
	void RejectStep();
	
public:
	StepControlSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp);
	
	FP Safety(const FP& v) const;

	static Vec<FP> GetTolerances(const Vec<FP>& a, const Vec<FP>& b, FP atol, FP rtol);
};

#endif
