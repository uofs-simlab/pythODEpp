#ifndef BACKWARD_EULER_H
#define BACKWARD_EULER_H

#include <core/common.h>
#include <methods/rk.h>

class BackwardEuler : public DIRK {
public:
	BackwardEuler(Hash<ParamValue>& params, BaseIVP* ivp) : DIRK(params, ivp, 1) {
		_a(0,0) = 1;
		_b(0) = 1;
		FillC();
	}

	virtual const char* GetName() const {
		return "Backward Euler";
	}
	
	virtual long GetOrder() const {
		return 1;
	}
};

#endif
