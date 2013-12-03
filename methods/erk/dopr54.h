#ifndef DOPR54_H
#define DOPR54_H

#include <core/common.h>
#include <methods/rk.h>

class DOPR54 : public ERK {
public:
	DOPR54(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 7) {
		_a(1,0) = 1./5;
		_a(2,0) = 3./40;  _a(2,1) = 9./40;
		_a(3,0) = 44./45; _a(3,1) = -56./15; _a(3,2) = 32./9;

		_a(4,0) = 19372./6561; _a(4,1) = -25360./2187; _a(4,2) = 64448./6561; _a(4,3) = -212./729;

		_a(5,0) = 9017./3168; _a(5,1) = -355./33;     _a(5,2) = 46732./5247;
		_a(5,3) = 49./176;    _a(5,4) = -5103./18656;

		_a(6,0) = 35./384;  _a(6,1) = 0;           _a(6,2) = 500./1113;
		_a(6,3) = 125./192; _a(6,4) = -2187./6784; _a(6,5) = 11./84;

		_b(0) = 35./384; _b(1) = 0; _b(2) = 500./1113; _b(3) = 125./192;
		_b(4) = -2187./6784; _b(5) = 11./84.;

		_baux(0) = 5179./57600; _baux(1) = 0; _baux(2) = 7571./16695; _baux(3) = 393./640;
		_baux(4) = -92097./339200; _baux(5) = 187./2100; _baux(6) = 1./40;
	
		FillC();
	}

	virtual const char* GetName() const {
		return "Dormand-Prince 5(4)";
	}
	
	virtual long GetOrder() const {
		return 5;
	}
	
	virtual long GetAuxOrder() const {
		return 4;
	}
};

#endif
