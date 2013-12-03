#ifndef NONSTIFF_E3_H
#define NONSTIFF_E3_H

class NonstiffE3 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = 0;
		jac(0,1) = 1;
		jac(1,0) = sqr(y[0])/2 - 1;
		jac(1,1) = 0;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = y[1];
		yp[1] = pow(y[0],3)/6 - y[0] + 2*sin(2.78535*t);
	}

public:
	NonstiffE3(Hash<ParamValue>& params) : BaseIVP(params) {
        params["tf"].SetFP(20);

        _initialCondition.Resize(2);
        _initialCondition[0] = 0;
        _initialCondition[1] = 0;
	}
		
	virtual const char* GetName() {
		return "Nonstiff E3";
	}
};

#endif

