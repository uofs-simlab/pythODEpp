#ifndef FORWARD_EULER_H
#define FORWARD_EULER_H

#include <core/common.h>
#include <methods/rk.h>

class ForwardEuler : public ERK {
public:
	ForwardEuler(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 1) {
		_b(0) = 1;
	}

	virtual const char* GetName() const {
		return "Forward Euler";
	}
	
	virtual long GetOrder() const {
		return 1;
	}
};

#endif
