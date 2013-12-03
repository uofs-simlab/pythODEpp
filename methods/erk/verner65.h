#ifndef VERNER65_H
#define VERNER65_H

#include <core/common.h>
#include <methods/rk.h>

class Verner65 : public ERK {
public:
	Verner65(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 8) {
		_a(1,0) = 1./6.;
		_a(2,0) = 4./75.;    _a(2,1) = 16./75.;
		_a(3,0) = 5./6.;     _a(3,1) = -8./3.; _a(3,2) = 5./2.;
		_a(4,0) = -165./64.; _a(4,1) = 55./6.; _a(4,2) = -425./64.;  _a(4,3) = 85./96.;
		_a(5,0) = 12./5.;    _a(5,1) = -8.;    _a(5,2) = 4015./612.; _a(5,3) = -11./36.; _a(5,4) = 88./255.;

		_a(6,0) = -8263./15000.; _a(6,1) = 124./75.;  _a(6,2) = -643./680.;     _a(6,3) = -81./250.;   _a(6,4) = 2484./10625.;
		_a(7,0) = 3501./1720.;   _a(7,1) = -300./43.; _a(7,2) = 297275./52632.; _a(7,3) = -319./2322.; _a(7,4) = 24068./84065.;
		_a(7,5) = 0.; _a(7,6) = 3850./26703.;



		_b(0) = 3./40.; _b(2) = 875./2244.;
		_b(3) = 23./72.; _b(4) = 264./1955.;
		_b(6) = 125./11592.; _b(7) = 43./616.;

		_baux(0) = 13./160.;
		_baux(2) = 2375./5984.;
		_baux(3) = 5./16.;
		_baux(4) = 12./85.;
		_baux(5) = 3./44.;

		FillC();
	}

	virtual const char* GetName() const {
		return "Verner 6(5)";
	}
	
	virtual long GetOrder() const {
		return 6;
	}
	
	virtual long GetAuxOrder() const {
		return 5;
	}
};

#endif
