#ifndef NONSTIFF_E4_H
#define NONSTIFF_E4_H

class NonstiffE4 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = 0;
		jac(0,1) = 1;
		jac(1,0) = 0;
		jac(1,1) = 0.8*y[1];
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = y[1];
		yp[1] = 0.032 - 0.4*sqr(y[1]);
	}

public:
	NonstiffE4(Hash<ParamValue>& params) : BaseIVP(params) {
        params["tf"].SetFP(20);

        _initialCondition.Resize(2);
        _initialCondition[0] = 30;
        _initialCondition[1] = 0;
	}
		
	virtual const char* GetName() {
		return "Nonstiff E4";
	}
};

#endif

