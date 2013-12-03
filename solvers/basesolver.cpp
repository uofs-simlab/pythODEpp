#include <solvers/basesolver.h>

void BaseSolver::CheckMaxSteps() {
	if( _steps >= _maxSteps )
		throw Exception() << "Maximum number of steps reached (" << _maxSteps << ").";

	if( _dt < 1e-10 )
		throw Exception() << "Minimum timestep reached.";
}

void BaseSolver::WriteFile(long f) {
	if( _lastWriteTime + _minWriteTime > _tn && !_complete && f )
		return;

	_lastWriteTime = _tn;
	
	std::ostringstream strs;
	strs << std::setfill('0') << std::setw(6) << f;
	std::string filename = _outPath + "/" + strs.str();
	
	std::ofstream file;
	file.open(filename.c_str(), std::ios::binary);
	
	if( !file.is_open() )
		throw Exception() << "Unable to open " << filename << ".";

	file.write((char*)&_tn, sizeof(_tn));    
	file << _yn;
	file.close();
}

void BaseSolver::UpdateTimestep() {
	_yn = _ynew;
	_tn += _dt;

	_method->UpdateTimestep();
}

BaseSolver::BaseSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp) : _ivp(ivp), _method(method) {
	if( !params.Get("path") )
		throw Exception() << "solver requires output path.";
	
	if( !params.Get("tf") )
		throw Exception() << "solver requires final time.";
	
	if( !params.Get("dt") )
		params["dt"].SetFP(1e-3);
	
	if( !params.Get("last step stretch") )
		params["last step stretch"].SetFP(1.1);

	if( !params.Get("max steps") )
		params["max steps"].SetLong(100000);

	if( _method )
		_method->SetSolverVariables(&_acceptedSteps, &_dtOld);

	_outPath = std::string(params["path"].GetString());
	_stretch = params["last step stretch"].GetFP();
	_maxSteps = params["max steps"].GetLong();
	
	_tn = params["tn"].GetFP();
	_dt = params["dt"].GetFP();
	_tf = params["tf"].GetFP();

	_printTime = (bool)params["print time"].GetLong();

	_complete = false;
	_lastWriteTime = _tn;
	_minWriteTime = 0;
	if( params.Get("min write time") )
		_minWriteTime = params["min write time"].GetFP();
		
	if( params.Get("timing group") && params["timing group"].GetFP() < 0 )
		throw Exception() << "timing group must be greater than or equal to zero.";

	if( _dt <= 0 )
		throw Exception() << "timestep cannot be zero.";
	
	_ivp->GetInitialCondition(_yn);
	_ynew.Resize(_yn.Size());

	_steps = 0;
	_acceptedSteps = 0;
	_rejectedSteps = 0;
	_dtOld = _dt;
}

BaseSolver::~BaseSolver() {
}

void BaseSolver::DumpRunInfo(Hash<ParamValue>& params) {
	// Set naming keys
	params["ivp name"].SetString(_ivp->GetName());
	params["method name"].SetString(_method->GetName());
	params["solver name"].SetString(GetName());

	params["time"].SetFP(_timer.msec());
	params["steps"].SetLong(_steps);
	params["accepted steps"].SetLong(_acceptedSteps);
	params["rejected steps"].SetLong(_rejectedSteps);

	_method->GetStats(params);
	_ivp->GetStats(params);

	std::string filename = _outPath + "/.runinfo";
	std::ofstream file;
	file.open(filename.c_str());

	if( !file.is_open() )
		throw Exception() << "Unable to write run information to file " << filename << ".";

	Hash<ParamValue>::Iterator it(params);
	while( *it ) {
		file << it.CurrentKey() << ":" << (**it).GetString() << std::endl;
		it.Next();
	}
	file.close();

	if( params.Get("print stats") && params["print stats"].GetLong() ) {
		std::cout << "Steps: " << _steps << "\n";
		std::cout << "Time:  " << (int)_timer.msec() << "ms\n";
	}
}

FP StepControlSolver::StandardStepControl() {
	long q = std::min(_method->GetOrder(), _method->GetAuxOrder());	
	return pow(_eps, FP(1)/(q+1));
}

FP StepControlSolver::PredictiveStepControl() {
	if( !_acceptedSteps || _rejectStep )
		return sqrt(_eps);
	
	if( !_epsLast )
		return 1/_maxChange;

	long q = std::min(_method->GetOrder(), _method->GetAuxOrder());	
	return (_dtOld/_dt)*pow(sqr(_eps)/_epsLast, FP(1)/(q+1));
}

FP StepControlSolver::CalculateFactor() {
    FP factor;
	switch( _scType ) {
	case STANDARD:
		factor = StandardStepControl();
		break;
	case PREDICTIVE:
		factor = PredictiveStepControl();
		break;
	case MIXED:
		factor = std::max(StandardStepControl(),
						  PredictiveStepControl());
		break;
	}

    return Safety(factor);
}

bool StepControlSolver::CheckMethodReject() {
	if( _ynew.IsNan() || !_method->Accept() ) {
		_rejectStep = true;
		_complete = false;
		_rejectedSteps++;
		_dt *= _minChange;
		return true;
	}

	return false;
}

void StepControlSolver::AcceptStep() {
	if( CheckMethodReject() )
		return;

	_rejectStep = false;
	UpdateTimestep();
	FP factor = CalculateFactor();
	_epsLast = _eps;
	_dtOld = _dt;
	_dt /= factor;
	_acceptedSteps++;
}

void StepControlSolver::RejectStep() {
	if( CheckMethodReject() )
		return;
	_rejectStep = true;
	_complete = false;
	_dt /= std::max(CalculateFactor(), 1/_maxRejectedChange);
	_rejectedSteps++;
}

StepControlSolver::StepControlSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp) : BaseSolver(params, method, ivp) {
	_rejectStep = false;
	_rTol = GetDefaultFP(params,"rtol",1e-5);
	_aTol = GetDefaultFP(params,"atol",1e-5);
	_safety = GetDefaultFP(params,"safety",0.9);
	_minChange = GetDefaultFP(params,"minchange",0.2);
	_maxChange = GetDefaultFP(params,"maxchange",5.);
	_maxRejectedChange = GetDefaultFP(params,"maxrejchange",1.);
	_restrictReject = (bool)GetDefaultLong(params,"restrict reject",1);
	
	ParamValue* pv = params.Get("step control");
	_scType = STANDARD;
	if( pv ) {
		if( std::string(pv->GetString()) == "Standard" )
			_scType = STANDARD;
		else if( std::string(pv->GetString()) == "Predictive" )
			_scType = PREDICTIVE;
		else if( std::string(pv->GetString()) == "Mixed" )
			_scType = MIXED;
		else
			throw Exception() << "Unknown step control type " << pv->GetString() << ".";
	}
}

FP StepControlSolver::Safety(const FP& v) const {
	if( !_restrictReject && _rejectStep )
		return v/_safety;

	return std::max(1/_maxChange, std::min(1/_minChange, v/_safety));
}

Vec<FP> StepControlSolver::GetTolerances(const Vec<FP>& a, const Vec<FP>& b, FP atol, FP rtol) {
	Vec<FP> ret(a.Size());
	for( long i = 0; i < a.Size(); i++ ) {
		FP ai = fabs(a[i]);
		FP bi = fabs(b[i]);
		ret[i] = (ai > bi ? ai : bi)*rtol+atol;
	}
	return ret;
}

