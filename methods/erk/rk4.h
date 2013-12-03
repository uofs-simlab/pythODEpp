#ifndef RK4_H
#define RK4_H

#include <core/common.h>
#include <methods/rk.h>

class RK4 : public ERK {
public:
	RK4(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 4) {
		_a(1,0) = 1./2;
		_a(2,1) = 1./2;
		_a(3,2) = 1.;
		_b(0) = 1./6;  _b(1) =1./3;
		_b(2) =1./3;   _b(3) = 1./6;
		FillC();
	}

	virtual const char* GetName() const {
		return "Runge-Kutta 4";
	}
	
	virtual long GetOrder() const {
		return 4;
	}
};

#endif
