#ifndef BS23_H
#define BS23_H

#include <core/common.h>
#include <methods/rk.h>

class BS23 : public ERK {
public:
	BS23(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 4) {
		_a(1,0) = 1./2;
		_a(2,0) = 0;  _a(2,1) = 3./4;
		_a(3,0) = 2./9; _a(3,1) = 1./3; _a(3,2) = 4./9;

		_b(0) = 2./9; _b(1) = 1./3; _b(2) = 4./9; _b(3) = 0;
		
		_baux(0) = 7./24; _baux(1) = 1./4; _baux(2) = 1./3; _baux(3) = 1./8;
	
		FillC();
	}

	virtual const char* GetName() const {
		return "Bogacki-Shampine 3(2)";
	}
	
	virtual long GetOrder() const {
		return 3;
	}
	
	virtual long GetAuxOrder() const {
		return 2;
	}
};

#endif

