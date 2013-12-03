#ifndef NONSTIFF_F3_H
#define NONSTIFF_F3_H

class NonstiffF3 : public BaseIVP {
protected:
	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = y[1];
		yp[1] = 0.01*y[1]*(1-sqr(y[0])) - y[0] - fabs(sin(M_PI*t));
	}

public:
	NonstiffF3(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(2);
		_initialCondition.Zero();
	}
		
	virtual const char* GetName() {
		return "Nonstiff F3";
	}
};

#endif

