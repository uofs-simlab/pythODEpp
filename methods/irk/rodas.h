#ifndef RADAS_H
#define RADAS_H

#include <core/common.h>
#include <core/csrmat.h>
#include <methods/basemethod.h>
#include <solvers/basesolver.h>

class RODAS : public BaseMethod {
	Mat<FP> _A, _C;
	Vec<FP> _c, _d;
	FP _gamma;

	Vec<FP>* _k;

public:
	RODAS(Hash<ParamValue>& params, BaseIVP* ivp) : BaseMethod(params, ivp) {
		_A.Resize(6,6);
		_C.Resize(6,6);
		_c.Resize(6);
		_d.Resize(6);
		
		_k = new Vec<FP>[6];
		if( ivp )
			for( long i = 0; i < 6; i++ )
				_k[i].Resize(ivp->Size());

		switch( params["elliptic"].GetLong() ) {
		case 0:
			_c(0) = 0;
			_c(1) = 0.386;
			_c(2) = 0.21;
			_c(3) = 0.63;
			_c(4) = 1;

			_d(0) =  0.2500000000000000e+0;
			_d(1) = -0.1043000000000000e+0;
			_d(2) =  0.1035000000000000e+0;
			_d(3) = -0.3620000000000023e-1;
			_d(4) =  0;

			_A(1,0) =  0.1544000000000000e+1;
			_A(2,0) =  0.9466785280815826e+0;
			_A(2,1) =  0.2557011698983284e+0;
			_A(3,0) =  0.3314825187068521e+1;
			_A(3,1) =  0.2896124015972201e+1;
			_A(3,2) =  0.9986419139977817e+0;
			_A(4,0) =  0.1221224509226641e+1;
			_A(4,1) =  0.6019134481288629e+1;
			_A(4,2) =  0.1253708332932087e+2;
			_A(4,3) = -0.6878860361058950e+0;

			_C(1,0) = -0.5668800000000000e+1;
			_C(2,0) = -0.2430093356833875e+1;
			_C(2,1) = -0.2063599157091915e+0;
			_C(3,0) = -0.1073529058151375e+0;
			_C(3,1) = -0.9594562251023355e+1;
			_C(3,2) = -0.2047028614809616e+2;
			_C(4,0) =  0.7496443313967647e+1;
			_C(4,1) = -0.1024680431464352e+2;
			_C(4,2) = -0.3399990352819905e+2;
			_C(4,3) =  0.1170890893206160e+2;
			_C(5,0) =  0.8083246795921522e+1;
			_C(5,1) = -0.7981132988064893e+1;
			_C(5,2) = -0.3152159432874371e+2;
			_C(5,3) =  0.1631930543123136e+2;
			_C(5,4) = -0.6058818238834054e+1;
			_gamma  =  0.2500000000000000e+0;
			break;

		case 1:
			_c(0) = 0;
			_c(1) = 0.3507221e0;
			_c(2) = 0.2557041e0;
			_c(3) = 0.6817790e0;
			_c(4) = 1;

			_d(0) =  0.2500000000000000e+0;
			_d(1) = -0.6902209999999998e-1;
			_d(2) = -0.9671999999999459e-3;
			_d(3) = -0.8797900000000025e-1;
			_d(4) =  0;

			_A(1,0) =  0.1402888400000000e+1;
			_A(2,0) =  0.6581212688557198e+0;
			_A(2,1) = -0.1320936088384301e+1;
			_A(3,0) =  0.7131197445744498e+1;
			_A(3,1) =  0.1602964143958207e+2;
			_A(3,2) = -0.5561572550509766e+1;
			_A(4,0) =  0.2273885722420363e+2;
			_A(4,1) =  0.6738147284535289e+2;
			_A(4,2) = -0.3121877493038560e+2;
			_A(4,3) =  0.7285641833203814e+0;

			_C(1,0) = -0.5104353600000000e+1;
			_C(2,0) = -0.2899967805418783e+1;
			_C(2,1) =  0.4040399359702244e+1;
			_C(3,0) = -0.3264449927841361e+2;
			_C(3,1) = -0.9935311008728094e+2;
			_C(3,2) =  0.4999119122405989e+2;
			_C(4,0) = -0.7646023087151691e+2;
			_C(4,1) = -0.2785942120829058e+3;
			_C(4,2) =  0.1539294840910643e+3;
			_C(4,3) =  0.1097101866258358e+2;
			_C(5,0) = -0.7629701586804983e+2;
			_C(5,1) = -0.2942795630511232e+3;
			_C(5,2) =  0.1620029695867566e+3;
			_C(5,3) =  0.2365166903095270e+2;
			_C(5,4) = -0.7652977706771382e+1;
			_gamma = 0.2500000000000000e+0; 
			break;

		case 2:
			_gamma = 0.25;
			_c(0) =  0;
			_c(1) =  3*_gamma;
			_c(2) =  0.21;
			_c(3) =  0.63;
			_c(4) =  1;

			_d(0) =  0.2500000000000000e+0;
			_d(1) = -0.5000000000000000e+0;
			_d(2) = -0.2350400000000000e-1;
			_d(3) = -0.3620000000000000e-1;

			_A(1,0) =  0.3000000000000000e+1;
			_A(2,0) =  0.1831036793486759e+1;
			_A(2,1) =  0.4955183967433795e+0;
			_A(3,0) =  0.2304376582692669e+1;
			_A(3,1) = -0.5249275245743001e-1;
			_A(3,2) = -0.1176798761832782e+1;
			_A(4,0) = -0.7170454962423024e+1;
			_A(4,1) = -0.4741636671481785e+1;
			_A(4,2) = -0.1631002631330971e+2;
			_A(4,3) = -0.1062004044111401e+1;

			_C(1,0) = -0.1200000000000000e+2;
			_C(2,0) = -0.8791795173947035e+1;
			_C(2,1) = -0.2207865586973518e+1;
			_C(3,0) =  0.1081793056857153e+2;
			_C(3,1) =  0.6780270611428266e+1;
			_C(3,2) =  0.1953485944642410e+2;
			_C(4,0) =  0.3419095006749676e+2;
			_C(4,1) =  0.1549671153725963e+2;
			_C(4,2) =  0.5474760875964130e+2;
			_C(4,3) =  0.1416005392148534e+2;
			_C(5,0) =  0.3462605830930532e+2;
			_C(5,1) =  0.1530084976114473e+2;
			_C(5,2) =  0.5699955578662667e+2;
			_C(5,3) =  0.1840807009793095e+2;
			_C(5,4) = -0.5714285714285717e+1;
			break;

		default:
			throw Exception() << "Rodas coefficients must be 0, 1, or 2.";
		}
	}

	virtual ~RODAS() {
		if( _k ) delete [] _k;
	}

	virtual void Step(FP tn, FP dt, const Vec<FP>& yn, Vec<FP>& ynew) {	
		BaseMat<FP>* dirmat;
		if( _sparse ) dirmat = new CSRMat<FP>(CSRMat<FP>::Eye(yn.Size())/(dt*_gamma) - *(CSRMat<FP>*)_ivp->JacSparse(tn,yn,1));
		else dirmat = new Mat<FP>(Mat<FP>::Eye(yn.Size())/(dt*_gamma) - *(Mat<FP>*)_ivp->Jac(tn,yn,1));
		dirmat->Factor();

		// Calculate df/dt for non-autonomous systems
		Vec<FP> dfdt(yn.Size());
		_ivp->RHSTimeDt(tn, yn, dfdt);

		// Stages 1 - 5
		Vec<FP> fn(yn.Size());
		for( long i = 0; i < 5; i++ ) {
			// Do function evaluation bit
			ynew = yn;
			for( long j = 0; j < i; j++ )
				ynew.AddScaled(_A(i,j),_k[j]);
			(*_ivp)(tn + dt*_c(i), ynew, fn);

			// Add non-autonomous term
			fn += dt*_d[i]*dfdt;

			// Add remaining stages
			for( long j = 0; j < i; j++ )
				fn.AddScaled(_C(i,j)/dt, _k[j]);

			// Solve for the new stage
			dirmat->Solve(fn, _k[i]);
		}
	
		// Add stage 5 (The coefficient on it must be 1,
		// and the rest of the coefficients must be the same as stage 4)
		ynew += _k[4];

		// Stage 6
		(*_ivp)(tn + dt, ynew, fn);
		for( long j = 0; j < 5; j++ )
			fn.AddScaled(_C(5,j)/dt, _k[j]);
		dirmat->Solve(fn, _k[5]);
		// Add stage 6 (The coefficient on it must be 1,
		// and the rest of the coefficients must be the same as stage 5)
		ynew += _k[5];

		delete dirmat;
	}

	virtual FP CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol) {
		return (_k[5]/StepControlSolver::GetTolerances(yn,ynew,atol,rtol)).RMS();
	}

	virtual const char* GetName() const {
		return "RODAS 4(3)";
	}
	
	virtual long GetOrder() const {
		return 4;
	}
	
	virtual long GetAuxOrder() const {
		return 3;
	}
};

#endif

