#ifndef NONSTIFF_F4_H
#define NONSTIFF_F4_H

class NonstiffF4 : public BaseIVP {
protected:
	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		if( t <= 10 )
			yp[0] = -2./21 - 120*(t-5)/(1 + 4*sqr(t-5));
		else
			yp[0] = -2*y[0];
	}

public:
	NonstiffF4(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(1);
		_initialCondition[0] = 1;
	}
		
	virtual const char* GetName() {
		return "Nonstiff F4";
	}
};

#endif

