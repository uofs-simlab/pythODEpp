#ifndef NONSTIFF_A4_H
#define NONSTIFF_A4_H

class NonstiffA4 : public BaseIVP {
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = (10-y(0))/40;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp(0) = (y(0)/4) * (1 - y(0)/20);
	}
	
public:
	NonstiffA4(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(1);
		_initialCondition[0] = 1;
	}
	
	virtual const char* GetName() {
		return "Nonstiff A4";
	}
};

#endif

