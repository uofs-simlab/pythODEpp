#ifndef NONSTIFF_F5_H
#define NONSTIFF_F5_H

class NonstiffF5 : public BaseIVP {
protected:
	FP _c;

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		FP pprime = 0;
		for( long i = 1; i < 20; i++ ) {
			CFP e = pow(CFP(t-i), 1./3);
			pprime += (t-i<0?-1:1)*(sqr(e.real())+sqr(e.imag()));
		}

		yp[0] = (4/(3*_c) * pprime * y[0]);
	}

public:
	NonstiffF5(Hash<ParamValue>& params) : BaseIVP(params), _c(0) {
		for( long i = 0; i < 20; i++ )
			_c += pow(i, 4./3);

		params["tf"].SetFP(20);

		_initialCondition.Resize(1);
		_initialCondition[0] = 1;
	}
		
	virtual const char* GetName() {
		return "Nonstiff F5";
	}
};

#endif

