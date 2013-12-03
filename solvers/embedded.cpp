#include <solvers/embedded.h>

EmbeddedSolver::EmbeddedSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp) : StepControlSolver(params, method, ivp) {
}

void EmbeddedSolver::RunSimulation() {
	WriteFile(0);

	_timer.Start();
	while( !_complete ) {
		do {
			CheckMaxSteps();
			
			_method->SetAccept(true);
			_method->PreStep(_tn, _dt, _yn);
			if( _tn + _dt*_stretch >= _tf ) {
				_complete = true;
				_dt = _tf - _tn;
			}

			if( _printTime )
				std::cout << "tn = " << _tn << " dt = " << _dt << std::endl;
		
			_method->Step(_tn, _dt, _yn, _ynew);
			_method->PostStep(_tn, _dt, _yn);
			_eps = _method->CalcEpsilon(_tn, _dt, _yn, _ynew, _aTol, _rTol);
			_eps <= 1 ? AcceptStep() : RejectStep();
			_steps++;
		} while( _rejectStep );

		WriteFile(_acceptedSteps);
	}
	_timer.End();

}

const char* EmbeddedSolver::GetName() const {
	return "Embedded Solver";
}

