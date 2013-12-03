#include <core/exception.h>
#include <core/csrmat.h>
#include <solvers/basesolver.h>
#include <methods/rk.h>

void RKMethod::FillC() {
	for( long i = 0; i < _m; i++ )
		for( long j = 0; j < _m; j++ )
			_c(i) += _a(i,j);
}

RKMethod::RKMethod(Hash<ParamValue>& params, BaseIVP* ivp, long m) : BaseMethod(params, ivp), _a(m,m), _b(m), _baux(m), _c(m), _m(m) {
	_k = new Vec<FP>[m];
	_a.Zero();
	_b.Zero();
	_c.Zero();
	_baux.Zero();
	
	// Allocate stages properly
	if( ivp )
		for( long i = 0; i < m; i++ )
			_k[i].Resize(ivp->Size());
}

const Mat<FP>& RKMethod::GetA() const {
	return _a;
}

const Vec<FP>& RKMethod::GetB() const {
	return _b;
}

const Vec<FP>& RKMethod::GetC() const {
	return _c;
}

FP RKMethod::CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol) {
	Vec<FP> aux(yn.Size());
	aux.Zero();
	for( long i = 0; i < _m; i++ )
		aux.AddScaled(_baux(i),_k[i]);
	aux *= dt;
	aux += yn;

	return ((ynew-aux)/StepControlSolver::GetTolerances(yn,ynew,atol,rtol)).RMS();
}

RKMethod::~RKMethod() {
	if( _k ) delete [] _k;
}

// -----------------------------------------------------------------------------------

ERK::ERK(Hash<ParamValue>& params, BaseIVP* ivp, long m) : RKMethod(params, ivp, m) {
}

void ERK::Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	ynew.Zero();

	for( long i = 0; i < _m; i++ ) {
		Vec<FP> k(yn.Size());
		k.Zero();
		
		for( long j = 0; j < i; j++ )
			k.AddScaled(_a(i,j), _k[j]);

		(*_ivp)(tn + dt*_c(i), yn + dt*k, _k[i]);
		ynew.AddScaled(_b(i), _k[i]);
	}
	
	ynew *= dt;
	ynew += yn;
}

// -----------------------------------------------------------------------------------

DIRK::DIRK(Hash<ParamValue>& params, BaseIVP* ivp, long m) : RKMethod(params, ivp, m) {
}

void DIRK::NewtonSolve(const FP tn, const FP dt, const long s, const Vec<FP>& yn, BaseMat<FP>* mat,
					   Vec<FP>& k, unsigned short split) {
	FP norm;	
	FP argt = tn + dt*_c(s);

	Vec<FP> f(yn.Size());
	Vec<FP> argy(yn.Size());

	for( long i = 0; i < 20; i++ ) {
		argy = yn;
		argy.AddScaled(dt*_a(s,s),k);
		(*_ivp)(argt, argy, f, split);
		f -= k;

		mat->Solve(f, argy);
		k += argy;

		norm = f.InfNorm();
		
		if( norm > _newtonFail )
			break;
		
		if( norm < _newtonTol || _ivp->JacobianSplitting() ) {
			return;
		}
	}

	// Method did not reach tolerence, so reject
	_accept = false;	
}

void DIRK::Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	ynew.Zero();
	
	Vec<FP> guess(yn.Size());
	(*_ivp)(tn, yn, guess);
	
	const BaseMat<FP>* jac = _sparse ? _ivp->JacSparse(tn,yn) : _ivp->Jac(tn,yn);
	_ivp->FreezeJacobian(true);
	
	BaseMat<FP>* mat = 0;
	
	for( long i = 0; i < _m; i++ ) {
		Vec<FP> accum(yn.Size());
		accum.Zero();
		
		for( long j = 0; j < i; j++ )
			accum += _a(i,j)*_k[j];
		accum *= dt;
		accum += yn;

		_k[i] = guess;
	
		if( _a(i,i) == 0 ) {
			(*_ivp)(tn + dt*_c(i), accum, _k[i]);
		} else if( i != 0 && _a(i,i) == _a(i-1, i-1) ) {
			NewtonSolve(tn, dt, i, accum, mat, _k[i]);
		} else {
			if( mat )
				delete mat;

			if( _sparse )
				mat = new CSRMat<FP>(CSRMat<FP>::Eye(yn.Size()) - dt*_a(i,i)**(CSRMat<FP>*)jac);
			else
				mat = new Mat<FP>(Mat<FP>::Eye(yn.Size()) - dt*_a(i,i)**(Mat<FP>*)jac);

			mat->Factor();
			NewtonSolve(tn, dt, i, accum, mat, _k[i]);
		}

		ynew.AddScaled(_b(i),_k[i]);
	}
	
	if( mat )
		delete mat;
	
	ynew *= dt;
	ynew += yn;
	_ivp->FreezeJacobian(false);
}

// -----------------------------------------------------------------------------------

void IMEX::FillC2() {
	for( long i = 0; i < _m; i++ )
		for( long j = 0; j < _m; j++ )
			_c2(i) += _a2(i,j);
}

IMEX::IMEX(Hash<ParamValue>& params, BaseIVP* ivp, long m) : DIRK(params, ivp, m), _a2(m,m), _b2(m), _baux2(m), _c2(m) {
	_k2 = new Vec<FP>[m];
	_a2.Zero();
	_b2.Zero();
	_baux2.Zero();
	_c2.Zero();
	
	// Allocate stages properly
	if( ivp )
		for( long i = 0; i < m; i++ )
			_k2[i].Resize(ivp->Size());
}

IMEX::~IMEX() {
	if( _k2 )
		delete [] _k2;
}

const Mat<FP>& IMEX::GetA2() const {
	return _a2;
}

const Vec<FP>& IMEX::GetB2() const {
	return _b2;
}

const Vec<FP>& IMEX::GetC2() const {
	return _c2;
}

void IMEX::Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	ynew.Zero();

	Timer jactimer;
	const BaseMat<FP>* jac = _sparse ? _ivp->JacSparse(tn,yn,1) : _ivp->Jac(tn,yn,1);
	if( _benchmark )
		printf("jac time: %dms\n", (int)jactimer.msec());
	_ivp->FreezeJacobian(true);
	
	BaseMat<FP>* mat = 0;
	Vec<FP> accum(yn.Size());
	Vec<FP> guess(yn.Size());
	
	(*_ivp)(tn, yn, guess, 1);

	Timer stagetimer;
	for( long i = 0; i < _m; i++ ) {
//		Timer st;
		accum.Zero();
		
		for( long j = 0; j < i; j++ ) {
			accum.AddScaled(_a(i,j),_k[j]);
			accum.AddScaled(_a2(i,j),_k2[j]);
		}

		accum *= dt;
		accum += yn;

		// Solve the implicit part
		_k[i] = guess;

		if( _a(i,i) == 0 ) {
			(*_ivp)(tn + dt*_c(i), accum, _k[i], 1);
		} else if( i != 0 && _a(i,i) == _a(i-1, i-1) ) {
			NewtonSolve(tn, dt, i, accum, mat, _k[i], 1);
		} else {
			if( mat )
				delete mat;

			if( _sparse ) {
				mat = new CSRMat<FP>(CSRMat<FP>::Eye(yn.Size()) - dt*_a(i,i)**(CSRMat<FP>*)jac);
			} else {
				mat = new Mat<FP>(Mat<FP>::Eye(yn.Size()) - dt*_a(i,i)**(Mat<FP>*)jac);
			}

		//	Timer ft;
			mat->Factor();
		//	printf("    factor time: %dms\n", (int)ft.msec());
			NewtonSolve(tn, dt, i, accum, mat, _k[i], 1);
		}

		accum.AddScaled(dt*_a(i,i), _k[i]);
		(*_ivp)(tn + dt*_c2(i), accum, _k2[i],2);

		// Add to the final vectors
		ynew.AddScaled(_b(i),_k[i]);
		ynew.AddScaled(_b2(i),_k2[i]);
		//printf(" stage %d: %dms\n", (int)i, (int)st.msec());
	}
	if( _benchmark )
		printf("stages: %dms\n", (int)stagetimer.msec());
	
	if( mat )
		delete mat;

	ynew *= dt;
	ynew += yn;
	_ivp->FreezeJacobian(false);
}

FP IMEX::CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol) {
	Vec<FP> auxF(yn.Size());
	Vec<FP> auxG(yn.Size());
	auxF.Zero();
	auxG.Zero();
	for( long i = 0; i < _m; i++ )
		auxF.AddScaled(_baux(i), _k[i]);
	for( long i = 0; i < _m; i++ )
		auxG.AddScaled(_baux2(i), _k2[i]);
	Vec<FP> aux = auxF+auxG;
	aux *= dt;
	aux += yn;
	return ((ynew-aux)/StepControlSolver::GetTolerances(yn,ynew,atol,rtol)).RMS();
}

