#ifndef NONSTIFF_A2_H
#define NONSTIFF_A2_H

class NonstiffA2 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = -3*y(0)*y(0)/2;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp(0) = -pow(y(0),3)/2;
	}
	
public:
	NonstiffA2(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(1);
		_initialCondition[0] = 1;
	}
	
	virtual const char* GetName() {
		return "Nonstiff A2";
	}
};

#endif

