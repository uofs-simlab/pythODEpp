#ifndef NONSTIFF_A1_H
#define NONSTIFF_A1_H

class NonstiffA1 : public BaseIVP {	
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = -1;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp(0) = -y(0);
	}

public:
	NonstiffA1(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(1);
		_initialCondition[0] = 1;
	}
	
	virtual const char* GetName() {
		return "Nonstiff A1";
	}
};

#endif

