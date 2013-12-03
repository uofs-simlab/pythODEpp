#ifndef NONSTIFF_F1_H
#define NONSTIFF_F1_H

class NonstiffF1 : public BaseIVP {
protected:
	FP _a;

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		FP pi2 = M_PI*M_PI;
		FP a2 = _a*_a;

		yp[0] = y[1];
		yp[1] = 2*_a*y[1] - (pi2+a2)*y[0] + (int(t)%2==0 ? 1 : -1);
	}

public:
	NonstiffF1(Hash<ParamValue>& params) : BaseIVP(params), _a(0.1) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(2);
		_initialCondition.Zero();
	}
		
	virtual const char* GetName() {
		return "Nonstiff F1";
	}
};

#endif

