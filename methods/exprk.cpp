#include <methods/exprk.h>

void AdditiveExpRK::ExpMtv(const BaseMat<FP>* M, const FP t, const Vec<FP>& v, Vec<FP>& expMtv) {
	Vec<FP> temp(v.Size());
	Vec<FP> tempMv = v;
	expMtv = v;

	for( long i = 1; i <= 3; i++ ) {
		M->VectorMult(tempMv,temp);
		tempMv = temp*(t/i);
		expMtv += tempMv;
	}
}

AdditiveExpRK::AdditiveExpRK(Hash<ParamValue>& params, BaseIVP* ivp) : BaseMethod(params, ivp) {
	bool flipExp = (bool)GetDefaultLong(params, "flipexp", 0);
	if( flipExp ) {
		_classical = 2;
		_exponential = 1;
	} else {
		_classical = 1;
		_exponential = 2;
	}
}

AdditiveExpRK::~AdditiveExpRK() {
}

// ------------------------------------------------------------------------------

DIRKCF1::DIRKCF1(Hash<ParamValue>& params, BaseIVP* ivp) : AdditiveExpRK(params, ivp) {
}

DIRKCF1::~DIRKCF1() {
}

void DIRKCF1::Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	if( _ivp->JacobianSplitting() )
		_ivp->Jac(tn, yn);
	_ivp->FreezeJacobian(true);

	const BaseMat<FP>* splitMat = _sparse ? _ivp->SplitMatSparse(tn, yn, _exponential) : _ivp->SplitMat(tn, yn, _exponential);

	Vec<FP> g(yn.Size());
	(*_ivp)(tn, yn, g, _classical);

	ExpMtv(splitMat, dt, yn, ynew);
	ynew += dt*g;

	_ivp->FreezeJacobian(false);
}

const char* DIRKCF1::GetName() const {
	return "DIRK-CF(1)";
}

long DIRKCF1::GetOrder() const {
	return 1;
}

// ------------------------------------------------------------------------------

DIRKCF2::DIRKCF2(Hash<ParamValue>& params, BaseIVP* ivp) : AdditiveExpRK(params, ivp) {
}

DIRKCF2::~DIRKCF2() {
}

void DIRKCF2::Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	// Set up Jacobian for implicit part
	BaseMat<FP>* jac;
	if( _sparse )
		jac = new CSRMat<FP>(*(CSRMat<FP>*)_ivp->JacSparse(tn, yn, _classical)*dt/2 - CSRMat<FP>::Eye(yn.Size()));
	else
		jac = new Mat<FP>(*(Mat<FP>*)_ivp->Jac(tn, yn, _classical)*dt/2 - Mat<FP>::Eye(yn.Size()));
	jac->Factor();
	_ivp->FreezeJacobian(true);

	// Explicit stage 1/2
	const BaseMat<FP>* expmat;
	if( _sparse )
		expmat = _ivp->SplitMatSparse(tn, yn, _exponential);
	else
		expmat = _ivp->SplitMat(tn, yn, _exponential);	
	Vec<FP> k_exp(yn.Size());
	ExpMtv(expmat, dt/2, yn, k_exp);

	// Implicit stage 1/2
	Vec<FP> k_imp = yn;
	Vec<FP> split1(yn.Size());
	for( int i = 0;; i++ ) {
		(*_ivp)(tn+dt/2, k_imp, split1, _classical);
		Vec<FP> f = k_exp + (dt/2)*split1 - k_imp;

		jac->Solve(f, split1);
		k_imp -= split1;
		
		FP norm = f.InfNorm();
		if( norm < _newtonTol )
			break;
		
		if( norm > _newtonFail || i > 10 ) {
			_accept = false;
			ynew.Zero();
			delete jac;
			return;
		}
	}

	// Complete method
	(*_ivp)(tn+dt/2, k_imp, split1, _classical);
	ExpMtv(expmat, -dt/2, split1, k_exp);

	if( _sparse )
		expmat = _ivp->SplitMatSparse(tn+dt, k_imp, _exponential);
	else
		expmat = _ivp->SplitMat(tn+dt, k_imp, _exponential);
	ExpMtv(expmat, dt, yn + dt*k_exp, ynew);

	_ivp->FreezeJacobian(false);
	delete jac;
}

const char* DIRKCF2::GetName() const {
	return "DIRK-CF(2)";
}

long DIRKCF2::GetOrder() const {
	return 2;
}

// ------------------------------------------------------------------------------

DIRKCF3::DIRKCF3(Hash<ParamValue>& params, BaseIVP* ivp) : AdditiveExpRK(params, ivp) {
}

DIRKCF3::~DIRKCF3() {
}

bool DIRKCF3::NewtonSolve(BaseMat<FP>* jac, Vec<FP>& guess, FP t, FP dt, const Vec<FP>& constant) {
	Vec<FP> temp(guess.Size());

	for( long i = 0;; i++ ) {
		(*_ivp)(t+dt, guess, temp, _classical);
		Vec<FP> f = constant + dt*temp - guess;
		jac->Solve(f, temp);
		guess -= temp;

		FP norm = f.InfNorm();
		if( norm < _newtonTol )
			break;

		if( norm > _newtonFail || i > 10 )
			return false;
	}

	return true;
}

void DIRKCF3::Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	FP gamma = (3.+sqrt(3.))/6;
	FP phi = 1./(6*(2*gamma-1));

	Vec<FP> constant2(yn.Size()), constant3(yn.Size());
	Vec<FP> y2g(yn.Size()), y3g(yn.Size());
	Vec<FP> p2ig2(yn.Size()), p3ig3(yn.Size());
	Vec<FP> y2(yn.Size()), y3(yn.Size());
	CSRMat<FP> p2m, p3m, p4m1, p4m2;

	const CSRMat<FP>* y1m;
	const CSRMat<FP>* y2m;
	const CSRMat<FP>* y3m;

	if( !_sparse )
		throw Exception() << "DIRKCF(3) does not support dense matrices.";

	yn.PrintMatlab(std::cout << "yn = ");
	// Set up Jacobian for implicit part
	CSRMat<FP>* jac = new CSRMat<FP>(*(CSRMat<FP>*)_ivp->JacSparse(tn, yn, _classical)*dt*gamma - CSRMat<FP>::Eye(yn.Size()));
	jac->Factor();
	_ivp->FreezeJacobian(true);

	// Calculate phi2 matrix
	std::cout << "exponential = " << _exponential << "\n";
	y1m = (const CSRMat<FP>*)_ivp->SplitMatSparse(tn, yn, _exponential);
	p2m = dt*gamma*(*y1m);
	p2m.ToDense().PrintMatlab(std::cout << "phi2 matrix = ");

	// Calculate stage 2
	ExpMtv(&p2m, 1, yn, constant2);
	if( !NewtonSolve(jac, y2 = yn, tn, gamma*dt, constant2) ) goto reject;
	y2.PrintMatlab(std::cout << "y2 = ");	

	// Calculate phi3 matrix
	y2m = (const CSRMat<FP>*)_ivp->SplitMatSparse(tn+gamma*dt, y2, _exponential);
	p3m = dt*((gamma-1)*(*y1m) + 2*(1-gamma)*(*y2m));
	p3m.ToDense().PrintMatlab(std::cout << "phi3 matrix = ");

	// Calculate stage 3
	(*_ivp)(tn + gamma*dt, y2, y2g, _classical);
	ExpMtv(&p2m, -1, y2, p2ig2);
	ExpMtv(&p3m, 1, yn + dt*(1-2*gamma)*p2ig2, constant3);
	constant3.PrintMatlab(std::cout << "constant3 = ");	
	if( !NewtonSolve(jac, y3 = y2, tn, gamma*dt, constant3) ) goto reject;
	y3.PrintMatlab(std::cout << "y3 = ");	
	
	// Calculate phi4 matrix
	y3m = (const CSRMat<FP>*)_ivp->SplitMatSparse(tn+(1-gamma)*dt, y3, _exponential);
	p4m1 = dt*(phi*(*y2m) - phi*(*y3m));
	p4m2 = dt*((1./2-phi)*(*y2m) + (1./2-phi)*(*y3m));

	// Calculate ynew
	(*_ivp)(tn + (1-gamma)*dt, y3, y3g, _classical);
	ExpMtv(&p3m, -1, y3g, p3ig3);
	ExpMtv(&p4m1, dt, yn + (dt/2)*(p2ig2 + p3ig3), ynew);
	ExpMtv(&p4m2, dt, ynew, ynew);

	_ivp->FreezeJacobian(false);
	delete jac;
	throw Exception() << "Done!\n";
	return;

reject:
	_accept = false;
	delete jac;
	ynew.Zero();
	throw Exception() << "Done!\n";
}

const char* DIRKCF3::GetName() const {
	return "DIRK-CF(3)";
}

long DIRKCF3::GetOrder() const {
	return 3;
}

// ------------------------------------------------------------------------------

ERKCF2::ERKCF2(Hash<ParamValue>& params, BaseIVP* ivp) : AdditiveExpRK(params, ivp) {
}

ERKCF2::~ERKCF2() {
}

void ERKCF2::Step(const FP tn, const FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	if( _ivp->JacobianSplitting() )
		_sparse ? _ivp->JacSparse(tn, yn) : _ivp->Jac(tn, yn);
	_ivp->FreezeJacobian(true);

	// Explicit stage 1/2
	const BaseMat<FP>* expmat;
	if( _sparse )
		expmat = _ivp->SplitMatSparse(tn, yn, _exponential);
	else
		expmat = _ivp->SplitMat(tn, yn, _exponential);
	Vec<FP> k_exp(yn.Size());
	ExpMtv(expmat, dt/2, yn, k_exp);

	// Finish first stage
	Vec<FP> split1(yn.Size());
	(*_ivp)(tn+dt/2, k_exp, split1, _classical);
	k_exp += (dt/2)*split1;

	// Complete method
	(*_ivp)(tn+dt/2, k_exp, split1, _classical);
	ExpMtv(expmat, -dt/2, split1, ynew);

	if( _sparse ) 
		expmat = _ivp->SplitMatSparse(tn+dt, k_exp, _exponential);
	else
		expmat = _ivp->SplitMat(tn+dt, k_exp, _exponential);
	ExpMtv(expmat, dt, yn + dt*ynew, ynew);

	_ivp->FreezeJacobian(false);
}

const char* ERKCF2::GetName() const {
	return "ERK-CF(2)";
}

long ERKCF2::GetOrder() const {
	return 2;
}

