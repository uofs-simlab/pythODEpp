#ifndef HEUN3_H
#define HEUN3_H

#include <core/common.h>
#include <methods/rk.h>

class Heun3 : public ERK {
public:
	Heun3(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 3) {
		_a(1,0) = 1./3; _a(2,1) = 2./3;
		_b(0) = 1./4;  _b(2) =3./4;
	
		FillC();
	}

	virtual const char* GetName() const {
		return "Heun 3";
	}
	
	virtual long GetOrder() const {
		return 3;
	}
};

#endif
