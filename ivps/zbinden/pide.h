#ifndef PIDE_H
#define PIDE_H

#include <ivps/splitivp.h>

class PIDE : public TwoSplittingIVP {
	long _n;
	FP _diff;
	FP _fac;
	Vec<FP> _y4;

public:
	PIDE(Hash<ParamValue>& params) : TwoSplittingIVP(params), _n(100), _diff(1.), _fac(1e-3) {
		params["tf"].SetFP(1);

		ParamValue* pv;
		if( (pv = params.Get("diff")) )
			_diff = pv->GetFP();	
		if( (pv = params.Get("fac")) )
			_fac = pv->GetFP();
		if( (pv = params.Get("N")) )
			_n = pv->GetFP();
			
		_initialCondition.Resize(_n);
		for( long i = 0; i < _n; i++ )
			_initialCondition[i] = 0.5*(1+cos(M_PI*(i+1)/_n));
		
		_y4.Resize(_n);
	}

	void Split1(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = _diff*(0.5*(2-sqrt(t)) - 2*y[0] + y[1]);
		for( long i = 1; i < _n-1; i++ )
			yp[i] = _diff*(y[i-1] - 2*y[i] + y[i+1]);
		yp[_n-1] = 2*_diff*(y[_n-2] - y[_n-1]);
	}

	void Split2(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		FP u04 = sqr(sqr( 0.5*(2-sqrt(t)) ));
		for( long i = 0; i < _n; i++ )
			_y4[i] = sqr(sqr(y[i]));

		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n-1; j++ )
				yp[i] -= _y4[j]/sqr(1+abs(FP(i+1)/_n - FP(j+1)/_n));
			yp[i] -= (u04/sqr(1+FP(i+1)/_n) + _y4[_n-1]/sqr(2-FP(i+1)/_n))/2;
		}

		yp *= _fac/_n;
	}

	virtual const char* GetName() {
		return "Parabolic integro-differential equation";
	}
};

#endif

