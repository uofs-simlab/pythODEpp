#include <solvers/constant.h>
#include <solvers/embedded.h>
#include <solvers/stepdoubling.h>

#define SOLVERCASE(solverclass) if( solver == #solverclass ) return new solverclass(params,method,ivp);

BaseIVP* AllocIVP(Hash<ParamValue>& params);
BaseMethod* AllocMethod(Hash<ParamValue>& params, BaseIVP* ivp);
BaseSolver* AllocSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp) {
	std::string solver = params["solver"].GetString();
	SOLVERCASE(ConstantSolver)
	SOLVERCASE(EmbeddedSolver)
	SOLVERCASE(StepDoublingSolver)
	throw Exception() << "Solver " << solver << " has not been defined.";
}

