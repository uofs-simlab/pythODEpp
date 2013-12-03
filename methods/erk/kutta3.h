#ifndef KUTTA3_H
#define KUTTA3_H

#include <core/common.h>
#include <methods/rk.h>

class Kutta3 : public ERK {
public:
	Kutta3(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 3) {
		_a(1,0) = 1./2;
		_a(2,0) = -1;
		_a(2,1) = 2;
		
		_b(0) = 1./6; _b(1) = 2./3; _b(2) = 1./6;
		FillC();
	}

	virtual const char* GetName() const {
		return "Kutta 3";
	}
	
	virtual long GetOrder() const {
		return 3;
	}
};

#endif
