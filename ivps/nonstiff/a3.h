#ifndef NONSTIFF_A3_H
#define NONSTIFF_A3_H

class NonstiffA3 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = cos(t);
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp(0) = y(0)*cos(t);
	}

public:
	NonstiffA3(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(1);
		_initialCondition[0] = 1;
	}
	
	virtual const char* GetName() {
		return "Nonstiff A3";
	}
};

#endif

