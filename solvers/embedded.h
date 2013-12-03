#ifndef EMBEDDED_SOLVER_H
#define EMBEDDED_SOLVER_H

#include <core/common.h>
#include <solvers/basesolver.h>

class EmbeddedSolver : public StepControlSolver {
public:
	EmbeddedSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp);
	
	virtual void RunSimulation();
	virtual const char* GetName() const;
};

#endif

