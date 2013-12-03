#ifndef RADAU5_H
#define RADAU5_H

#include <core/common.h>
#include <methods/basemethod.h>
#include <solvers/basesolver.h>

class Radau5 : public BaseMethod {
	Mat<FP> _a;
	Mat<FP> _Tr;
	Mat<FP> _Ti;
	Vec<FP> _c;
	Vec<FP> _d;

	Vec<FP> _cont1;
	Vec<FP> _cont2;
	Vec<FP> _cont3;
	Vec<FP> _Z1;
	Vec<FP> _Z2;
	Vec<FP> _Z3;

	BaseMat<FP>* _E1;
	BaseMat<CFP>* _E2;
	
	FP _gamma;
	FP _alpha;
	FP _beta;

	FP _newtonFail;
	FP _newtonTol;
	
public:
	Radau5(Hash<ParamValue>& params, BaseIVP* ivp) : BaseMethod(params, ivp), _a(3,3), _Tr(3,3), _Ti(3,3), _c(2), _d(3), _newtonFail(1e20), _newtonTol(1e-8) {
		FP sq6 = sqrt(6);
		
		_a(0,0) = (88-7*sq6)/360;
		_a(0,1) = (296-169*sq6)/1800;
		_a(0,2) = (-2+3*sq6)/225;

		_a(1,0) = (296+169*sq6)/1800;
		_a(1,1) = (88+7*sq6)/360;
		_a(1,2) = (-2-3*sq6)/225;

		_a(2,0) = (16-6*sq6)/36;
		_a(2,1) = (16+6*sq6)/36;
		_a(2,2) = 1./9;
	
		_c(0) = (4-sq6)/10;
		_c(1) = (4+sq6)/10;

		_d(0) = (-13-7*sq6)/3;
		_d(1) = (-13+7*sq6)/3;
		_d(2) = -1./3;

		_gamma = 30/(6 + pow(81,1./3) - pow(9,1./3));
		_alpha = (12 - pow(81,1./3) + pow(9,1./3))/60;
		_beta  = (pow(81,1./3) + pow(9,1./3))*sqrt(3)/60;

		FP mag = sqr(_alpha) + sqr(_beta);
		_alpha /= mag;
		_beta  /= mag;

		_Tr(0,0) = 9.1232394870892942792e-2;
		_Tr(0,1) = -0.14125529502095420843;
		_Tr(0,2) = -3.0029194105147424492e-2;
		_Tr(1,0) = 0.24171793270710701896;
		_Tr(1,1) = 0.20412935229379993199;
		_Tr(1,2) = 0.38294211275726193779;
		_Tr(2,0) = 0.96604818261509293619;
		_Tr(2,1) = 1;
		_Tr(2,2) = 0;

		_Ti(0,0) = 4.3255798900631553510;
		_Ti(0,1) = 0.33919925181580986954;
		_Ti(0,2) = 0.54177053993587487119;
		_Ti(1,0) = -4.1787185915519047273;
		_Ti(1,1) = -0.32768282076106238708;
		_Ti(1,2) = 0.47662355450055045196;
		_Ti(2,0) = -0.50287263494578687595;
		_Ti(2,1) = 2.5719269498556054292;
		_Ti(2,2) = -0.59603920482822492497;

		ParamValue* pv;
		if( (pv = params.Get("newton fail")) )
			_newtonFail = pv->GetFP();
		if( (pv = params.Get("newton tol")) )
			_newtonTol = pv->GetFP();

		if( _ivp ) {
			_Z1.Resize(_ivp->Size());
			_Z2.Resize(_ivp->Size());
			_Z3.Resize(_ivp->Size());
		}
		
		if( _sparse ) {
			_E1 = new CSRMat<FP>;
			_E2 = new CSRMat<CFP>;
		} else {
			_E1 = new Mat<FP>;
			_E2 = new Mat<CFP>;
		}
	}

	~Radau5() {
		delete _E1;
		delete _E2;
	}

	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {
		FP g = _gamma/dt;
		FP a = _alpha/dt;
		FP b = _beta/dt;

		if( _sparse ) {
			CSRMat<FP>* jac = (CSRMat<FP>*)_ivp->JacSparse(tn,yn);
			*(CSRMat<FP>*)_E1 = g*CSRMat<FP>::Eye(yn.Size()) - *jac;
			*(CSRMat<CFP>*)_E2 = CFP(a,b)*CSRMat<CFP>::Eye(yn.Size()) - CSRMat<CFP>(*jac);
		} else {
			Mat<FP>* jac = (Mat<FP>*)_ivp->Jac(tn,yn);
			*(Mat<FP>*)_E1 = g*Mat<FP>::Eye(yn.Size()) - *jac;
			*(Mat<CFP>*)_E2 = CFP(a,b)*Mat<CFP>::Eye(yn.Size()) - *jac;
		}

		_E1->Factor();
		_E2->Factor();		

		Vec<FP> F1(yn.Size());
		Vec<FP> F2(yn.Size());
		Vec<FP> F3(yn.Size());

		if( !*_acceptedSteps ) {
			_Z1.Zero(); _Z2.Zero(); _Z3.Zero();
			F1.Zero(); F2.Zero(); F3.Zero();
		} else {
			FP c3q = dt / *_dtOld;
			FP c1q = _c(0)*c3q;
			FP c2q = _c(1)*c3q;

			FP c2m1 = _c(1) - 1;
			FP c1m1 = _c(0) - 1;
			Vec<FP> Z1I = c1q*(_cont1+(c1q-c2m1)*(_cont2+(c1q-c1m1)*_cont3));
			Vec<FP> Z2I = c2q*(_cont1+(c2q-c2m1)*(_cont2+(c2q-c1m1)*_cont3));
			Vec<FP> Z3I = c3q*(_cont1+(c3q-c2m1)*(_cont2+(c3q-c1m1)*_cont3));
			_Z1 = Z1I;
			_Z2 = Z2I;
			_Z3 = Z3I;
			F1 = _Ti(0,0)*Z1I + _Ti(0,1)*Z2I + _Ti(0,2)*Z3I;
			F2 = _Ti(1,0)*Z1I + _Ti(1,1)*Z2I + _Ti(1,2)*Z3I;
			F3 = _Ti(2,0)*Z1I + _Ti(2,1)*Z2I + _Ti(2,2)*Z3I;
		}

		for( long i = 0; i < 20; i++ ) {
			Vec<FP> A1(yn.Size()), A2(yn.Size()), A3(yn.Size());

			(*_ivp)(tn + _c(0)*dt, yn + _Z1, A1);
			(*_ivp)(tn + _c(1)*dt, yn + _Z2, A2);
			(*_ivp)(tn + dt, yn + _Z3, A3);
			
			_Z1 = _Ti(0,0)*A1 + _Ti(0,1)*A2 + _Ti(0,2)*A3;
			_Z2 = _Ti(1,0)*A1 + _Ti(1,1)*A2 + _Ti(1,2)*A3;
			_Z3 = _Ti(2,0)*A1 + _Ti(2,1)*A2 + _Ti(2,2)*A3;	

			_Z1 -= F1*g;
			_Z2 += -F2*a + F3*b;
			_Z3 += -F3*a - F2*b;

			Vec<FP> temp(_Z1);
			_E1->Solve(temp, _Z1);

			Vec<CFP> res(yn.Size());
			Vec<CFP> tempC(Vec<CFP>(_Z2) + CFP(0,1)*Vec<CFP>(_Z3));
			_E2->Solve(tempC, res);
			
			_Z2 = VecReal(res);
			_Z3 = VecImag(res);

			FP norm = 0;
			for( long j = 0; j < yn.Size(); j++ )
				norm += sqr(_Z1(j)) + sqr(_Z2(j)) + sqr(_Z3(j));
			norm = sqrt(norm/(3*yn.Size()));

			// Do the newton update
			F1 += _Z1;
			F2 += _Z2;
			F3 += _Z3;
			_Z1 = _Tr(0,0)*F1 + _Tr(0,1)*F2 + _Tr(0,2)*F3;
			_Z2 = _Tr(1,0)*F1 + _Tr(1,1)*F2 + _Tr(1,2)*F3;
			_Z3 = _Tr(2,0)*F1 + _Tr(2,1)*F2 + _Tr(2,2)*F3;

			if( norm > _newtonFail )
				break;

			if( norm < _newtonTol ) {
				ynew = yn + _Z3;
				return;
			}
		}

		_accept = false;
	}

	virtual FP CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol) {
		Vec<FP> fn(yn.Size());
		(*_ivp)(tn, yn, fn);
		Vec<FP> diff = (_d(0)*_Z1 + _d(1)*_Z2 +_d(2)*_Z3)/dt;

		// First prediction
		Vec<FP> err(yn.Size());
		Vec<FP> temp(fn+diff);
		_E1->Solve(temp, err);
		err /= StepControlSolver::GetTolerances(yn,yn,atol,rtol);
		FP eps = err.RMS();
	
		//if( eps < 1 && *_acceptedSteps )
		if( eps < 1 )
			return eps;

		// Second prediction
		(*_ivp)(tn, err + yn, fn);
		temp = fn + diff;
		_E1->Solve(temp, err);
		err /= StepControlSolver::GetTolerances(yn,yn,atol,rtol);
	
		return err.RMS();
	}

	virtual void UpdateTimestep() {
		BaseMethod::UpdateTimestep();
		
		Vec<FP> ak = (_Z1-_Z2)/(_c(0)-_c(1));
		_cont1 = (_Z2 - _Z3)/(_c(1)-1);
		_cont2 = (ak - _cont1)/(_c(0)-1);
		_cont3 = _cont2 - (ak-_Z1/_c(0))/_c(1);
	}

	virtual const char* GetName() const {
		return "Radau IIA 5";
	}
	
	virtual long GetOrder() const {
		return 5;
	}

	virtual long GetAuxOrder() const {
		return 3;
	}
};

#endif

