#ifndef RUNGE3_H
#define RUNGE3_H

#include <core/common.h>
#include <methods/rk.h>

class Runge3 : public ERK {
public:
	Runge3(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 4) {
		_a(1,0) = 1./2;
		_a(2,1) = 1;
		_a(3,2) = 1;
		
		_b(0) = 1./6; _b(1) = 2./3; _b(3) = 1./6;
		FillC();
	}

	virtual const char* GetName() const {
		return "Runge 3";
	}
	
	virtual long GetOrder() const {
		return 3;
	}
};

#endif
