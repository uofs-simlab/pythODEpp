#ifndef MERSON43_H
#define MERSON43_H

#include <core/common.h>
#include <methods/rk.h>

class Merson43 : public ERK {
public:
	Merson43(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 5) {
		_a(1,0) = 1./3;
		_a(2,0) = 1./6; _a(2,1) = 1./6;
		_a(3,0) = 1./8; _a(3,1) = 0;    _a(3,2) =  3./8;
		_a(4,0) = 1./2; _a(4,1) = 0;    _a(4,2) = -3./2; _a(4,3) = 2.;

		_b(0) = 1./6; _b(1) = 0; _b(2) = 0; _b(3) = 2./3; _b(4) = 1./6;

		_baux(0) = 1./10; _baux(1) = 0; _baux(2) = 3./10; _baux(3) = 2./5; _baux(4) = 1./5;
	
		FillC();
	}

	virtual const char* GetName() const {
		return "Merson 4(3)";
	}
	
	virtual long GetOrder() const {
		return 4;
	}
	
	virtual long GetAuxOrder() const {
		return 3;
	}
};

#endif
