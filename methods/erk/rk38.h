#ifndef RK38_H
#define RK38_H

#include <core/common.h>
#include <methods/rk.h>

class RK38 : public ERK {
public:
	RK38(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 4) {
		_a(1,0) =  1./3;
		_a(2,0) = -1./3; _a(2,1) =  1.;
		_a(3,0) =  1.;   _a(3,1) = -1.; _a(3,2) = 1.;

		_b(0) = 1./8;  _b(1) = 3./8;
		_b(2) = 3./8;  _b(3) = 1./8;
		FillC();
	}

	virtual const char* GetName() const {
		return "Three-Eighths Rule 4";
	}
	
	virtual long GetOrder() const {
		return 4;
	}
};

#endif

