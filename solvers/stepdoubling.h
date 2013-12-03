#ifndef STEP_DOUBLING_SOLVER_H
#define STEP_DOUBLING_SOLVER_H

#include <core/common.h>
#include <solvers/basesolver.h>

class StepDoublingSolver : public StepControlSolver {
	Vec<FP> _ystep2;

	virtual void UpdateTimestep();

public:
	StepDoublingSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp);
	
	virtual void RunSimulation();
	virtual const char* GetName() const;
};

#endif

