#ifndef RUNGE2_H
#define RUNGE2_H

#include <core/common.h>
#include <methods/rk.h>

class Runge2 : public ERK {
public:
	Runge2(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 2) {
		_a(1,0) = 1./2;
		_b(1) = 1.;
		FillC();
	}

	virtual const char* GetName() const {
		return "Runge 2";
	}
	
	virtual long GetOrder() const {
		return 2;
	}
};

#endif
