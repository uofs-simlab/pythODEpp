#ifndef CONSTANT_SOLVER_H
#define CONSTANT_SOLVER_H

#include <core/common.h>
#include <solvers/basesolver.h>

class ConstantSolver : public BaseSolver {
public:
	ConstantSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp);
	
	virtual void RunSimulation();
	virtual const char* GetName() const;
};

#endif
