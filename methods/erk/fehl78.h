#ifndef FEHL78_H
#define FEHL78_H

#include <core/common.h>
#include <methods/rk.h>

class FEHL78 : public ERK {
public:
	FEHL78(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 13) {
		_a(1,0)  = 2./27.;
		_a(2,0)  = 1./36.; _a(2,1) = 1./12.;
		_a(3,0)  = 1./24.; _a(3,2) = 1./8.;
		_a(4,0)  = 5./12.; _a(4,2) = -25./16.; _a(4,3) = 25./16.;

		_a(5,0)  = 1./20.; _a(5,3) = 1./4.; _a(5,4) = 1./5.;

		_a(6,0)  = -25./108.; _a(6,3) = 125./108.; _a(6,4) = -65./27.; _a(6,5) = 125./54.;

		_a(7,0)  = 31./300.; _a(7,4) = 61./225.;
		_a(7,5) = -2./9.; _a(7,6) = 13./900.;

		_a(8,0)  = 2.;       _a(8,3) = -53./6.; _a(8,4) = 704./45.;
		_a(8,5)  = -107./9.; _a(8,6) = 67./90.; _a(8,7) = 3.; 

		_a(9,0)  = -91./108.; _a(9,3) = 23./108.; _a(9,4) = -976./135.; _a(9,5) = 311./54.; _a(9,6) = -19./60.; _a(9,7) = 17./6.; _a(9,8) = -1./12.;

		_a(10,0) = 2383./4100.; _a(10,3) = -341./164.;  _a(10,4) = 4496./1025.;
		_a(10,5) = -301./82.;   _a(10,6) = 2133./4100.; _a(10,7) = 45./82.;
		_a(10,8) = 45./164.;    _a(10,9) = 18./41.;

		_a(11,0) = 3./205.; _a(11,5) = -6./41.; _a(11,6) = -3./205.;
		_a(11,7) = -3./41.; _a(11,8) = 3./41.;  _a(11,9) = 6./41.;

		_a(12,0) = -1777./4100.; _a(12,3) = -341./164.;  _a(12,4) = 4496./1025.;
		_a(12,5) = -289./82.;    _a(12,6) = 2193./4100.; _a(12,7) = 51./82.;
		_a(12,8) = 33./164.;     _a(12,9) = 12./41.;     _a(12,11) = 1.; 

		_b(0) = 41./840.; _b(5) = 34./105.; _b(6) = 9./35.;
		_b(7) = 9./35.; _b(8) = 9./280.; _b(9) = 9./280.;
		_b(10) = 41./840.;

		_baux(5) = 34./105.; _baux(6) = 9./35.;
		_baux(7) = 9./35.;
		_baux(8) = 9./280.;
		_baux(9) = 9./280.;
		_baux(11) = 41./840.;
		_baux(12) = 41./840.;
		
		FillC();
	}

	virtual const char* GetName() const {
		return "Runge-Kutta-Fehlberg 7(8)";
	}
	
	virtual long GetOrder() const {
		return 7;
	}
	
	virtual long GetAuxOrder() const {
		return 8;
	}
};

#endif

