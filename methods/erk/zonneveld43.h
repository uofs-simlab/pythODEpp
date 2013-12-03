#ifndef ZONNEVELD43_H
#define ZONNEVELD43_H

#include <core/common.h>
#include <methods/rk.h>

class Zonneveld43 : public ERK {
public:
	Zonneveld43(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 5) {
		_a(1,0) = 1./2;  _a(2,1) = 1./2;  _a(3,2) = 1.;
		_a(4,0) = 5./32; _a(4,1) = 7./32; _a(4,2) = 13./32; _a(4,3) = -1./32;

		_b(0) = 1./6; _b(1) = 1./3; _b(2) = 1./3; _b(3) = 1./6;

		_baux(0) = -1./2; _baux(1) = 7./3; _baux(2) = 7./3;
		_baux(3) = 13./6; _baux(4) = -16./3;
	
		FillC();
	}

	virtual const char* GetName() const {
		return "Zonneveld 4(3)";
	}
	
	virtual long GetOrder() const {
		return 4;
	}
	
	virtual long GetAuxOrder() const {
		return 3;
	}
};

#endif
