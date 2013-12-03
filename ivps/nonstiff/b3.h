#ifndef NONSTIFF_B3_H
#define NONSTIFF_B3_H

class NonstiffB3 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = -1; jac(0,1) =  0;      jac(0,2) = 0;
		jac(1,0) =  1; jac(1,1) = -2*y[1]; jac(1,2) = 0;
		jac(2,0) =  0; jac(2,1) =  2*y[1]; jac(2,2) = 0;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = -y[0];
		yp[1] = y[0] - y[1]*y[1];
		yp[2] = y[1]*y[1];
	}

public:
	NonstiffB3(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(3);
		_initialCondition[0] = 1;
		_initialCondition[1] = 0;
		_initialCondition[2] = 0;
	}
	
	virtual const char* GetName() {
		return "Nonstiff B3";
	}
};

#endif

