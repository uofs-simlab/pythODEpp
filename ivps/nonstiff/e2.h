#ifndef NONSTIFF_E2_H
#define NONSTIFF_E2_H

class NonstiffE2 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = 0;
		jac(0,1) = 1;
		jac(1,0) = -2*y[0]*y[1] - 1;
		jac(0,1) = 1-sqr(y[0]);
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = y[1];
		yp[1] = (1-sqr(y[0]))*y[1] - y[0];
	}

public:
	NonstiffE2(Hash<ParamValue>& params) : BaseIVP(params) {
        params["tf"].SetFP(20);

        _initialCondition.Resize(2);
        _initialCondition[0] = 2;
        _initialCondition[1] = 0;
	}
		
	virtual const char* GetName() {
		return "Nonstiff E2";
	}
};

#endif

