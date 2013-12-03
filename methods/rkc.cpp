#include <core/exception.h>
#include <solvers/basesolver.h>
#include <methods/rkc.h>

FP RKC2::EstimateJacobianSpectralRadius(FP t, const Vec<FP>& y, Vec<FP>& guess, const Vec<FP>& yp, long split) {	
	if( !guess.Size() )
		guess = yp;
	
	FP yNorm = y.Norm();
	FP ypNorm = yp.Norm();
	FP eps = std::numeric_limits<FP>().epsilon();
	FP sqrtEps = sqrt(eps);

	FP dyNorm;
	if( yNorm != 0 && ypNorm != 0 ) {
		dyNorm = yNorm*sqrtEps;
		guess = y + guess*dyNorm/ypNorm;
	} else if( yNorm != 0 ) {
		dyNorm = yNorm*sqrtEps;
		guess = y + y*sqrtEps;
	} else if( ypNorm != 0 ) {
		dyNorm = eps;
		guess = guess*dyNorm/ypNorm;
	} else {
		dyNorm = eps;
		guess = Vec<FP>::Ones(guess.Size())*dyNorm;
	}

	// Iterate with nonlinear power method
	FP sigma = 0;
	FP spRad = 0;
	
	for( long i = 0; i < 50; i++ ) {
		Vec<FP> fv(y.Size());
		(*_ivp)(t, guess, fv, split);
				
		FP fvypNorm = (fv-yp).Norm();
		FP sigma1 = sigma;
		sigma = fvypNorm / dyNorm;
		spRad = 1.2*sigma;
		if( i > 1 && fabs(sigma - sigma1) < sigma*0.01 ) {
			guess -= y;
			break;
		}

		if( fvypNorm != 0 )
			guess = y + (fv-yp) * dyNorm / fvypNorm;
		else {
			long index = i % y.Size();
			guess[index] = y[index] - (guess[index] - y[index]);
		}
	}

	return spRad;
}

RKC2::RKC2(Hash<ParamValue>& params, BaseIVP* ivp) : BaseMethod(params, ivp) {
	if( params.Get("eta") )
		_eta = params["eta"].GetFP();
	else
		_eta = 2./13;

	if( params.Get("max stages") )
		_maxStages = params["max stages"].GetFP();
	else
		_maxStages = 100;

	_statMaxStages = 2;
	_statMinStages = _maxStages;
	_statSteps = 0;
	_statStages = 0;

	if( ivp )
		_F0.Resize(ivp->Size());
}

void RKC2::PreStep(const FP tn, FP& dt, Vec<FP>& yn) {
	// Estimate spectral radius of the Jacobian
	(*_ivp)(tn, yn, _F0);
	FP spRad = EstimateJacobianSpectralRadius(tn, yn, _spRadGuess, _F0);

	// Calculate the number of stages
	_m = 1 + (long)sqrt(dt*spRad/0.65 + 1);

	if( _m > _maxStages ) {
		_m = _maxStages;
		dt = 0.65*(_m*_m-1)/spRad;
	}
	
	if( _m < _statMinStages )
		_statMinStages = _m;
	if( _m > _statMaxStages )
		_statMaxStages = _m;
	_statStages += _m;
	_statSteps++;

}

void RKC2::Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	const Vec<FP>& K0 = yn;
	Vec<FP>& Km = ynew;

	Vec<FP> Kjm1(yn.Size());
	Vec<FP> Kjm2(yn.Size());
	Vec<FP> Fjm1(yn.Size());

	// Begin steps
	FP w0 = 1 + _eta/(_m*_m);
	FP w1 = Cheb1p(_m, w0) / Cheb1pp(_m, w0);

	// Start setting up constants
	FP bj = Cheb1ppRecursive(2, 2*w0*w0-1, 2*w0, w0) / sqr(Cheb1pRecursive(2,2*w0));
	FP bjm1 = bj;
	FP bjm2 = bj;

	// Calculate K1
	FP kj   = bjm1*w1;
	Kjm1 = K0;
	Km   = K0 + kj*dt*_F0;
	
	// Set up Chebyshev polynomials
	FP Tjm1 = 1;    // T0(w0) = 1
	FP Tj   = w0;   // T1(w0) = w0
	FP Ujm1 = 1;    // U0(w0) = 1
	FP Uj   = 2*w0; // U1(w0) = w0
	FP Tjm2, Ujm2;

	for( long j = 2; j < _m+1; j++ ) {
		// Update chebyshev polynomials
		Tjm2 = Tjm1;
		Tjm1 = Tj;
		Tj = 2*w0*Tjm1-Tjm2;
		Ujm2 = Ujm1;
		Ujm1 = Uj;
		Uj = 2*w0*Ujm1-Ujm2;

		// Update past stages
		Kjm2 = Kjm1;
		Kjm1 = Km;

		// Set up constants
		bjm2 = bjm1;
		bjm1 = bj;
		bj = Cheb1ppRecursive(j,Tj,Ujm1,w0) / sqr(Cheb1pRecursive(j,Ujm1));
		FP ajm1 = 1 - bjm1*Tjm1;
		FP muj = 2*bj*w0/bjm1;
		FP nuj = -bj/bjm2;
		kj = 2*bj*w1/bjm1;

		// Update function evaluations
		FP cjm1;
		if( j == 2 )
			cjm1 = w1*bjm1;
		else
			cjm1 = w1*Cheb1ppRecursive(j-1,Tjm1,Ujm2,w0)/Cheb1pRecursive(j-1,Ujm2);
		(*_ivp)(tn+cjm1*dt, Kjm1, Fjm1);
		
		// Update stage
		Km = (1-muj-nuj)*K0;
		Km.AddScaled(muj,Kjm1);
		Km.AddScaled(nuj,Kjm2);
		Km.AddScaled(kj*dt,Fjm1);
		Km.AddScaled(-ajm1*kj*dt,_F0);
	}
}

FP RKC2::CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol) {
	Vec<FP> f1(yn.Size());
	(*_ivp)(tn, yn, f1);
	
	Vec<FP> f2(yn.Size());
	(*_ivp)(tn+dt, ynew, f2);

	Vec<FP> error = 0.8*(yn-ynew) + 0.4*dt*(f1 + f2);
	error /= StepControlSolver::GetTolerances(yn, ynew, atol, rtol);
	return error.RMS();
}

void RKC2::GetStats(Hash<ParamValue>& params) const {
	BaseMethod::GetStats(params);
	params["max stages"].SetFP(_statMaxStages);
	params["min stages"].SetFP(_statMinStages);
	params["avg stages"].SetFP(FP(_statStages)/_statSteps);
}

const char* RKC2::GetName() const {
	return "Runge-Kutta-Chebyshev (2)";
}

long RKC2::GetOrder() const {
	return 2;
}

long RKC2::GetAuxOrder() const {
	return 2;
}

FP RKC2::Cheb1(long n, FP x) {
	if( x >= 1 )
		return cosh(n*acosh(x));
	if( x <= -1 )
		return ((n&1)?-1:1) * cosh(n*acosh(-x));
	return cos(n*acos(x));
}

FP RKC2::Cheb2(long n, FP x) {
	if( x ==  0 ) return Cheb1(n,x);
	if( x ==  1 ) return n + 1;
	if( x == -1 ) return ((n&1)?-1:1) * (n+1);
	if( x  >  1 ) return sinh((n+1)*acosh(x))/sinh(acosh(x));
	if( x  < -1 ) return ((n&1)?-1:1) * sinh((n+1)*acosh(-x))/sinh(acosh(-x));
	return sin((n+1)*acos(x))/sin(acos(x));
}

FP RKC2::Cheb1p(long n, FP x) {
	return n*Cheb2(n-1, x);
}

FP RKC2::Cheb1pp(long n, FP x) {
	if( x == 1 ) return (n*n*(1-n*n))/3.;
	if( x == -1 ) return ((n&1)?-1:1) * Cheb1pp(n,1);
	return (n*Cheb1(n,x) - x*Cheb2(n-1,x))*n/(x*x-1);
}

FP RKC2::Cheb1pRecursive(long n, FP Ujm1) {
	return n*Ujm1;
}

FP RKC2::Cheb1ppRecursive(long n, FP Tj, FP Ujm1, FP x) {
	if( x == 1 )
		return (n*n*(1-n*n))/3.;
	return (n*Tj - x*Ujm1)*n/(x*x-1);
}

// ----------------------------------------------------------------------------

RKC1::RKC1(Hash<ParamValue>& params, BaseIVP* ivp) : RKC2(params, ivp) {
}

void RKC1::Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	const Vec<FP>& K0 = yn;
	Vec<FP>& Km = ynew;

	Vec<FP> Kjm1(yn.Size());
	Vec<FP> Kjm2(yn.Size());
	Vec<FP> Fjm1(yn.Size());

	// Begin steps
	FP w0 = 1 + _eta/(_m*_m);
	FP w1 = Cheb1(_m, w0) / Cheb1p(_m, w0);
	
	// Set up Chebyshev polynomials
	FP Tjm1 = 1;    // T0(w0) = 1
	FP Tj   = w0;   // T1(w0) = w0
	FP Ujm1 = 1;    // U0(w0) = 1
	FP Uj   = 2*w0; // U1(w0) = w0
	FP Tjm2, Ujm2;

	// Calculate K1
	Kjm1 = K0;
	Km   = K0 + w1/w0*dt*_F0;
	
	for( long j = 2; j < _m+1; j++ ) {
		// Update chebyshev polynomials
		Tjm2 = Tjm1;
		Tjm1 = Tj;
		Tj = 2*w0*Tjm1-Tjm2;
		Ujm2 = Ujm1;
		Ujm1 = Uj;
		Uj = 2*w0*Ujm1-Ujm2;

		// Update past stages
		Kjm2 = Kjm1;
		Kjm1 = Km;

		// Set up constants
		FP bjm2 = 1/Tjm2;
		FP bjm1 = 1/Tjm1;
		FP bj   = 1/Tj;
		FP muj = 2*bj*w0/bjm1;
		FP nuj = -bj/bjm2;
		FP kj = 2*bj*w1/bjm1;

		// Update function evaluations
		FP cjm1 = w1*Cheb1pRecursive(j-1,Ujm2)/Tjm1;
		(*_ivp)(tn+cjm1*dt, Kjm1, Fjm1);

		// Update stage
		Km.Zero();
		Km.AddScaled(muj,Kjm1);
		Km.AddScaled(nuj,Kjm2);
		Km.AddScaled(kj*dt,Fjm1);
	}
}

const char* RKC1::GetName() const {
	return "Runge-Kutta-Chebyshev (1)";
}

// ----------------------------------------------------------------------------

PRKC::PRKC(Hash<ParamValue>& params, BaseIVP* ivp) : RKC2(params, ivp) {
	if( ivp ) {
		_Gm1.Resize(ivp->Size());
		_G0.Resize(ivp->Size());
		_Gmm1.Resize(ivp->Size());
		_Gm.Resize(ivp->Size());
		_K0.Resize(ivp->Size());
		_Km.Resize(ivp->Size());
		_Kf.Resize(ivp->Size());
	}
}

void PRKC::PreStep(const FP tn, FP& dt, Vec<FP>& yn) {
	// For jacobian splitting
	if( _ivp->JacobianSplitting() ) {
		_ivp->Jac(tn,yn);
		_ivp->FreezeJacobian(true);
	}

	// Estimate spectral radius of the Jacobian
	(*_ivp)(tn, yn, _F0, 1);
	FP spRadF = EstimateJacobianSpectralRadius(tn, yn, _spRadGuess, _F0, 1);

	(*_ivp)(tn, yn, _Gm1, 2);	
	FP spRadG = _ivp->JacobianSplitting() ? 0 : EstimateJacobianSpectralRadius(tn, yn, _spRadGuessG, _Gm1, 2);
	if( spRadG*dt > 1.7 ) dt = 1.7/spRadG;

	// Calculate the number of stages
	_m = 1 + (long)sqrt(dt*spRadF/0.65 + 1);
	
	if( _m > _maxStages ) {
		_m = _maxStages;
		dt = 0.65*(_m*_m-1)/spRadF;
	}

	if( _m < _statMinStages )
		_statMinStages = _m;
	if( _m > _statMaxStages )
		_statMaxStages = _m;
	_statStages += _m;
	_statSteps++;
	
	if( _ivp->JacobianSplitting() )
		_ivp->FreezeJacobian(false);
}

void PRKC::Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	// For jacobian splitting
	if( _ivp->JacobianSplitting() )
		_ivp->FreezeJacobian(true);

	// Calculate K0 right away
	FP alpha0 = 1./2;
	_K0 = yn + alpha0*dt*_Gm1;

	// Memory to store function evaluations and stages
	Vec<FP> Kjm1(yn.Size());
	Vec<FP> Kjm2(yn.Size());
	Vec<FP> Fjm1(yn.Size());

	// Begin steps
	FP w0 = 1 + _eta/(_m*_m);
	FP w1 = Cheb1p(_m, w0) / Cheb1pp(_m, w0);

	// Start setting up constants
	FP bj = Cheb1ppRecursive(2, 2*w0*w0-1, 2*w0, w0) / sqr(Cheb1pRecursive(2,2*w0));
	FP bjm1 = bj;
	FP bjm2 = bj;

	// Recalculate F0 from prestep because it isn't the same for PRKC
	(*_ivp)(tn, _K0, _F0, 1);

	// Calculate K1
	FP kj   = bjm1*w1;
	Kjm1 = _K0;
	_Kf  = _K0 + kj*dt*_F0;

	// Set up Chebyshev polynomials
	FP Tjm1 = 1;    // T0(w0) = 1
	FP Tj   = w0;   // T1(w0) = w0
	FP Ujm1 = 1;    // U0(w0) = 1
	FP Uj   = 2*w0; // U1(w0) = w0
	FP Tjm2, Ujm2;

	for( long j = 2; j < _m+1; j++ ) {
		// update chebyshev polynomials
		Tjm2 = Tjm1;
		Tjm1 = Tj;
		Tj = 2*w0*Tjm1-Tjm2;
		Ujm2 = Ujm1;
		Ujm1 = Uj;
		Uj = 2*w0*Ujm1-Ujm2;
		
		// Update past stages
		Kjm2 = Kjm1;
		Kjm1 = _Kf;
		
		// Set up constants
		bjm2 = bjm1;
		bjm1 = bj;
		bj = Cheb1ppRecursive(j,Tj,Ujm1,w0) / sqr(Cheb1pRecursive(j,Ujm1));
		FP ajm1 = 1 - bjm1*Tjm1;
		FP muj = 2*bj*w0/bjm1;
		FP nuj = -bj/bjm2;
		kj = 2*bj*w1/bjm1;

		// Update function evaluations
		if( j == 2 )
			_cmm1 = w1*bjm1;
		else
			_cmm1 = w1*Cheb1ppRecursive(j-1,Tjm1,Ujm2,w0)/Cheb1pRecursive(j-1,Ujm2);
		(*_ivp)(tn+_cmm1*dt, Kjm1, Fjm1, 1);
	
		// Update stage
		_Kf = (1-muj-nuj)*_K0;
		_Kf.AddScaled(muj,Kjm1);
		_Kf.AddScaled(nuj,Kjm2);
		_Kf.AddScaled(kj*dt,Fjm1);
		_Kf.AddScaled(-ajm1*kj*dt,_F0);
	}

	// mth stage
	FP alpha1 = -1.5;
	FP alpha2 = 2;
	FP alpha3 = 0;

	(*_ivp)(tn+alpha0*dt, _K0, _G0, 2);
	(*_ivp)(tn+alpha0*dt, Kjm1, _Gmm1, 2);
	
	_Km = _Kf;
	_Km.AddScaled(alpha1*dt, _Gm1);
	_Km.AddScaled(alpha2*dt, _G0);
	_Km.AddScaled(alpha3*dt, _Gmm1);
	
	// m+1th stage
	FP alpha4 = -1./3;
	FP alpha5 = (2.*_cmm1 - 1.) / (3.*_cmm1);
	FP alpha6 = 1./(3.*_cmm1);
	FP alpha7 = 1./6;
	(*_ivp)(tn+dt, _Km, _Gm, 2);

	ynew = _Kf;
	ynew.AddScaled(alpha4*dt, _Gm1);
	ynew.AddScaled(alpha5*dt, _G0);
	ynew.AddScaled(alpha6*dt, _Gmm1);
	ynew.AddScaled(alpha7*dt, _Gm);

	if( _ivp->JacobianSplitting() )
		_ivp->FreezeJacobian(false);
}

FP PRKC::CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol) {
	if( _ivp->JacobianSplitting() )
		_ivp->FreezeJacobian(true);
	
	Vec<FP> f1(yn.Size());
	(*_ivp)(tn, _K0, f1, 1);
	
	Vec<FP> f2(yn.Size());
	(*_ivp)(tn+dt, _Kf, f2, 1);
	
	if( _ivp->JacobianSplitting() )
		_ivp->FreezeJacobian(false);

	Vec<FP> errF = 0.8*(_K0-_Kf) + 0.4*dt*(f1 + f2);

	FP beta1 = -1./2;
	FP beta3 = 1/(2*_cmm1);
	FP beta2 = 1 - beta3;
	Vec<FP> errG = ynew - (_Kf + dt*(beta1*_Gm1 + beta2*_G0 + beta3*_Gmm1));
	
	errF /= StepControlSolver::GetTolerances(_K0, _Kf, atol, rtol);
	errG /= StepControlSolver::GetTolerances(yn, ynew, atol, rtol);
	return std::max(errF.RMS(), errG.RMS());
}

const char* PRKC::GetName() const {
	return "Partitioned Runge-Kutta-Chebyshev";
}

// ----------------------------------------------------------------------------

IRKC::IRKC(Hash<ParamValue>& params, BaseIVP* ivp) : RKC2(params, ivp) {
}

void IRKC::NewtonSolve(const Mat<FP>& LU, const Mat<FP>& P, const Vec<FP>& constant, FP t, FP dt, FP k1, Vec<FP>& k, Vec<FP>& Gj) {
    Vec<FP> f(k.Size());
    for( long i = 0; i < 20; i++ ) {
        (*_ivp)(t, k, Gj, 2);
		f  = -k1*dt*Gj;
		f -= constant;
        f += k;
        
		k -= Mat<FP>::SolveLU(LU,P,f);

        FP norm = f.InfNorm();
        if( norm > _newtonFail )
            break;
        else if( norm < _newtonTol )
            return;
    }       

    // Method did not reach tolerence, so reject
    _accept = false;    
}

void IRKC::PreStep(const FP tn, FP& dt, Vec<FP>& yn) {
	// Estimate spectral radius of the Jacobian
	(*_ivp)(tn, yn, _F0, 1);
	FP spRad = EstimateJacobianSpectralRadius(tn, yn, _spRadGuess, _F0, 1);

	// Calculate the number of stages
	_m = 1 + (long)sqrt(dt*spRad/0.65 + 1);
	_m = 2;
	
	if( _m > _maxStages ) {
		_m = _maxStages;
		dt = 0.65*(_m*_m-1)/spRad;
	}

	if( _m < _statMinStages )
		_statMinStages = _m;
	if( _m > _statMaxStages )
		_statMaxStages = _m;
	_statStages += _m;
	_statSteps++;
}

void IRKC::Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
	const Vec<FP>& K0 = yn;
	Vec<FP>& Km = ynew;

	Vec<FP> Kjm1(yn.Size());
	Vec<FP> Kjm2(yn.Size());
	Vec<FP> Fjm1(yn.Size());
	Vec<FP> Gj(yn.Size());
	Vec<FP> Gjm1(yn.Size());
	Vec<FP> Gjm2(yn.Size());

	// Begin steps
	FP w0 = 1 + _eta/(_m*_m);
	FP temp1 = w0*w0 - 1;
	FP temp2 = sqrt(temp1);
	FP arg = _m*log(w0 + temp2);
	FP w1 = sinh(arg) * temp1 / (cosh(arg)*_m*temp2 - w0*sinh(arg)); 

	// Start setting up constants
	FP bj   = 1/w0;
	FP bjm1 = 1/sqr(2*w0);
	FP bjm2;
	FP cj   = w1/w0;
	FP cjm1 = 0;
	FP cjm2 = 0;
	FP aj   = 0;
	FP ajm1 = 0;
	FP kj   = bj*w1;
	_k1 = kj;

	// Calculate Jacobian at t0
	_jac = (Mat<FP>*)_ivp->Jac(tn, yn, 2);
    Mat<FP> LU, P;
    (Mat<FP>::Eye(yn.Size()) - (_k1*dt)**_jac).CalcLU(LU,P);

	// Calculate FO and G0
	Vec<FP> G0(yn.Size());
	(*_ivp)(tn, yn, G0, 2);
	Gjm1 = G0;

	// Calculate K1
	Kjm1 = K0;
	Km   = K0; // Guess
	NewtonSolve(LU, P, K0 + _k1*dt*_F0, tn + cj*dt, dt, _k1, Km, Gj);

	// Set up Chebyshev polynomials
	FP Tjm1 = 1;    // T0(w0) = 1
	FP Tj   = w0;   // T1(w0) = w0
	FP Ujm1 = 1;    // U0(w0) = 1
	FP Uj   = 2*w0; // U1(w0) = w0
	FP Tjm2, Ujm2;

	for( long j = 2; j < _m+1; j++ ) {
		// Update chebyshev polynomials
		Tjm2 = Tjm1;
		Tjm1 = Tj;
		Tj = 2*w0*Tjm1-Tjm2;
		Ujm2 = Ujm1;
		Ujm1 = Uj;
		Uj = 2*w0*Ujm1-Ujm2;

		// Update past stages and G evaluations
		Kjm2 = Kjm1;
		Kjm1 = Km;
		Gjm2 = Gjm1;
		Gjm1 = Gj;

		// Set up constants
		bjm2 = bjm1; bjm1 = bj;
		bj = Cheb1ppRecursive(j,Tj,Ujm1,w0) / sqr(Cheb1pRecursive(j,Ujm1));
		FP muj = 2*bj*w0/bjm1;
		FP nuj = -bj/bjm2;
		kj = 2*bj*w1/bjm1;
		ajm1 = aj;
		aj = 1 - bj*Tj;
		cjm2 = cjm1; cjm1 = cj;
		cj = muj*cjm1 + nuj*cjm2 + aj*(1+kj);

		// Update function evaluations
		(*_ivp)(tn+cjm1*dt, Kjm1, Fjm1, 1);
		// Update stage
		Vec<FP> stage = (1-muj-nuj)*K0;
		stage.AddScaled(muj,Kjm1);
		stage.AddScaled(nuj,Kjm2);
		stage.AddScaled(kj*dt,Fjm1);
		stage.AddScaled(-ajm1*kj*dt,_F0);
		stage.AddScaled(-(ajm1*kj+(1-muj-nuj)*_k1)*dt,G0);
		stage.AddScaled(-nuj*_k1*dt,Gjm2);
		NewtonSolve(LU, P, stage, tn + cj*dt, dt, _k1, Km, Gj);
	}
}

FP IRKC::CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol) {
	Mat<FP> LU, P;
    (Mat<FP>::Eye(yn.Size()) - dt**_jac).CalcLU(LU,P);

	Vec<FP> fn1(yn.Size()), fn2(yn.Size());
	(*_ivp)(tn+dt, ynew, fn1);
	(*_ivp)(tn, yn, fn2);

	Vec<FP> rhs = (dt/2)*(fn1 - fn2);

	(*_ivp)(tn+dt, ynew, fn1, 2);
	(*_ivp)(tn, yn, fn2, 2);
	rhs += dt*_k1*(fn1 - fn2);

	Vec<FP> err = Mat<FP>::SolveLU(LU, P, rhs); 
	err /= StepControlSolver::GetTolerances(yn, ynew, atol, rtol);
	return err.RMS();
}


const char* IRKC::GetName() const {
	return "IMEX Runge-Kutta-Chebyshev";
}

