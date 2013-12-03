#ifndef RKF45_H
#define RKF45_H

#include <core/common.h>
#include <methods/rk.h>

class RKF45 : public ERK {
public:
	RKF45(Hash<ParamValue>& params, BaseIVP* ivp) : ERK(params, ivp, 6) {
		_a(1,0) = 1./4;
		_a(2,0) = 3./32;       _a(2,1) = 9./32;
		_a(3,0) = 1932./2197.; _a(3,1) = -7200./2197.; _a(3,2) = 7296./2197;
		_a(4,0) = 439./216.;   _a(4,1) = -8.;          _a(4,2) = 3680./513.;   _a(4,3) = -845./4104.;
		_a(5,0) = -8./27.;     _a(5,1) = 2.;           _a(5,2) = -3544./2565.; _a(5,3) = 1859./4104.; _a(5,4) = -11./40.;
		
		_b(0) = 25./216.; _b(1) = 0.; _b(2) = 1408./2565.; _b(3) = 2197./4104.; _b(4) = -1./5.; _b(5) = 0;
		
		_baux(0) = 16./135.; _baux(1) = 0.; _baux(2) = 6656./12825.;
		_baux(3) = 28561./56430.; _baux(4) = -9./50.; _baux(5) = 2./55;
		
		FillC();
	}

	virtual const char* GetName() const {
		return "Runge-Kutta-Fehlberg 4(5)";
	}
	
	virtual long GetOrder() const {
		return 4;
	}
	
	virtual long GetAuxOrder() const {
		return 5;
	}
};

#endif
