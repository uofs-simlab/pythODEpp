#ifndef NONSTIFF_E5_H
#define NONSTIFF_E5_H

class NonstiffE5 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = 0;
		jac(0,1) = 1;
		jac(1,0) = 0;
		jac(1,1) = y[1] / ( sqrt(1+sqr(y[1])) * (25 - t) );
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = y[1];
		yp[1] = sqrt(1+sqr(y[1]))/(25 - t);
	}

public:
	NonstiffE5(Hash<ParamValue>& params) : BaseIVP(params) {
        params["tf"].SetFP(20);

        _initialCondition.Resize(2);
        _initialCondition[0] = 0;
        _initialCondition[1] = 0;
	}
		
	virtual const char* GetName() {
		return "Nonstiff E5";
	}
};

#endif

