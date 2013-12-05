#ifndef ARK1_H
#define ARK1_H

#include <methods/rk.h>

class ARK1 : public IMEX {
public:
	ARK1(Hash<ParamValue>& params, BaseIVP* ivp) : IMEX(params, ivp, 1) {
		_a(0,0) = 1.;
		_b(0) = _b2(0) = 1.;

		FillC();
		FillC2();
	}

	virtual const char* GetName() const {
		return "ARK 1";
	}
	
	virtual long GetOrder() const {
		return 1;
	}
};

#endif

