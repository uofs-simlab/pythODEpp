#include <solvers/stepdoubling.h>

StepDoublingSolver::StepDoublingSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp) : StepControlSolver(params, method, ivp) {
	_ystep2.Resize(_yn.Size());
}

void StepDoublingSolver::RunSimulation() {
	WriteFile(0);

	_timer.Start();
	while( !_complete ) {
		do {
			CheckMaxSteps();
			
			_method->SetAccept(true);
			_method->PreStep(_tn, _dt, _yn);
			if( _tn + _dt*(1+_stretch) >= _tf ) {
				_complete = true;
				_dt = (_tf - _tn)/2;
			}

			if( _printTime )
				std::cout << "tn = " << _tn << " dt = " << 2*_dt << std::endl;

			// Do all the stepping
			_method->Step(_tn, _dt, _yn, _ystep2);
			_method->Step(_tn+_dt, _dt, _ystep2, _ynew);
			_method->Step(_tn, 2*_dt, _yn, _ystep2);
			_method->PostStep(_tn, _dt, _yn);

			Vec<FP> err = (_ynew - _ystep2)/(pow(2., (FP)_method->GetOrder())-1);
			_eps = (err/StepControlSolver::GetTolerances(_yn,_ynew,_aTol,_rTol)).RMS();
			_eps <= 1 ? AcceptStep() : RejectStep();
			_steps++;
		} while( _rejectStep );

		WriteFile(_acceptedSteps);
	}
	_timer.End();

}

void StepDoublingSolver::UpdateTimestep() {
    _yn = _ynew;
    _tn += 2*_dt;

    _method->UpdateTimestep();
}

const char* StepDoublingSolver::GetName() const {
	return "Step Doubling Solver";
}

