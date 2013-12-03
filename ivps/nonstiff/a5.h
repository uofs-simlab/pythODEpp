#ifndef NONSTIFF_A5_H
#define NONSTIFF_A5_H

class NonstiffA5 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = 2*t/((y(0)+t)*(y(0)+t));
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp(0) = (y(0)-t) / (y(0)+t);
	}
	
public:
	NonstiffA5(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(1);
		_initialCondition[0] = 4;
	}
	
	virtual const char* GetName() {
		return "Nonstiff A5";
	}
};

#endif

