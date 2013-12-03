#ifndef NONSTIFF_B2_H
#define NONSTIFF_B2_H

class NonstiffB2 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = -1; jac(0,1) =  1; jac(0,2) = 0;
		jac(1,0) =  1; jac(1,1) = -2; jac(1,2) = 1;
		jac(2,0) =  0; jac(2,1) =  1; jac(2,2) = -1;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = -y[0] + y[1];
		yp[1] = y[0] - 2*y[1] + y[2];
		yp[2] = y[1] - y[2];
	}
	
public:
	NonstiffB2(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(3);
		_initialCondition[0] = 2;
		_initialCondition[1] = 0;
		_initialCondition[2] = 1;
	}
	
	virtual const char* GetName() {
		return "Nonstiff B2";
	}
};

#endif

