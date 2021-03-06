#ifndef BS54_H
#define BS54_H

#include <core/common.h>
#include <methods/rk.h>

class BS54 : public ERK {
public:
	BS54(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 8) {
		_a(1,0) = 1./6;
		_a(2,0) = 2./27; _a(2,1) = 4./27;
		_a(3,0) = 183./1372; _a(3,1) = -162./343; _a(3,2) = 1053./1372;
		_a(4,0) = 68./297; _a(4,1) = -4./11; _a(4,2) = 42./143; _a(4,3) = 1960./3861;
		_a(5,0) = 597./22528; _a(5,1) = 81./352; _a(5,2) = 63099./585728;
			_a(5,3) = 58653./366080; _a(5,4) = 4617./20480;
		_a(6,0) = 174197./959244; _a(6,1) = -30942./79937; _a(6,2) = 8152137./19744439;
			_a(6,3) = 666106./1039181; _a(6,4) = -29421./29068; _a(6,5) = 482048./414219;
		_a(7,0) = 587./8064; _a(7,1) = 0; _a(7,2) = 4440339./15491840;
			_a(7,3) = 24353./124800; _a(7,4) = 387./44800; _a(7,5) = 2152./5985;
			_a(7,6) = 7267./94080;

		_b(0) = 587./8064; _b(1) = 0; _b(2) = 4440339./15491840; _b(3) = 24353./124800;
			_b(4) = 387./44800; _b(5) = 2152./5985; _b(6) = 7267./94080; _b(7) = 0;
		
		_baux(0) = 2479./34992; _baux(1) = 0; _baux(2) = 123./416;
			_baux(3) = 612941./3411720; _baux(4) = 43./1440; _baux(5) = 2272./6561;
			_baux(6) = 79937./1113912; _baux(7) = 3293./556956;
	
		FillC();
	}

	virtual const char* GetName() const {
		return "Bogacki-Shampine 5(4)";
	}
	
	virtual long GetOrder() const {
		return 5;
	}
	
	virtual long GetAuxOrder() const {
		return 4;
	}
};

#endif

