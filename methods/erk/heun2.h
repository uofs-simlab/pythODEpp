#ifndef HEUN2_H
#define HEUN2_H

#include <core/common.h>
#include <methods/rk.h>

class Heun2 : public ERK {
public:
	Heun2(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 2) {
		_a(1,0) = 1.;
		_b(0) = 1./2;  _b(1) =1./2;
		_baux(0) = 0.65; _baux(1) = 0.35;
	
		FillC();
	}

	virtual const char* GetName() const {
		return "Heun 2";
	}
	
	virtual long GetOrder() const {
		return 2;
	}
};

#endif
