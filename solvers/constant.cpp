#include <solvers/constant.h>

ConstantSolver::ConstantSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp) : BaseSolver(params, method, ivp) {
}

void ConstantSolver::RunSimulation() {
	WriteFile(0);

	_timer.Start();
	while( !_complete ) {
		CheckMaxSteps();
		_method->PreStep(_tn, _dt, _yn);
		
		if( _tn + _dt*_stretch >= _tf ) {
			_complete = true;
			_dt = _tf - _tn;
		}

		if( _printTime )
			printf("tn = %g, dt = %g\n", _tn, _dt);

		_method->Step(_tn, _dt, _yn, _ynew);
		_method->PostStep(_tn, _dt, _yn);

		UpdateTimestep();
		WriteFile(++_steps);
	}
	_timer.End();
}

const char* ConstantSolver::GetName() const {
	return "Constant Solver";
}

