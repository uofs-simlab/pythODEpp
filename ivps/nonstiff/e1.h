#ifndef NONSTIFF_E1_H
#define NONSTIFF_E1_H

class NonstiffE1 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		FP tp1 = t+1;

		jac(0,0) = 0;
		jac(0,1) = 1;
		jac(1,0) = -(1-0.25/sqr(tp1));
		jac(1,1) = -1/tp1;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		FP tp1 = t+1;

		yp[0] = y[1];
		yp[1] = -(y[1]/tp1 + (1-0.25/sqr(tp1))*y[0]); 
	}

public:
	NonstiffE1(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(2);
		_initialCondition[0] = 0.6713967071418031;
		_initialCondition[1] = 0.09540051444747454;
	}
		
	virtual const char* GetName() {
		return "Nonstiff E1";
	}
};

#endif

