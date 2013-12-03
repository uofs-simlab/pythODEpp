#ifndef NONSTIFF_B5_H
#define NONSTIFF_B5_H

class NonstiffB5 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = 0;          jac(0,1) = y[2];       jac(0,2) = y[1];
		jac(1,0) = -y[2];      jac(1,1) = 0;          jac(1,2) = -y[0];
		jac(2,0) = -0.51*y[1]; jac(2,1) = -0.51*y[0]; jac(2,2) = 0;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = y[1]*y[2];
		yp[1] = -y[0]*y[2];
		yp[2] = -0.51*y[0]*y[1];
	}

public:
	NonstiffB5(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(3);
		_initialCondition[0] = 0;
		_initialCondition[1] = 1;
		_initialCondition[2] = 1;
	}
	
	virtual const char* GetName() {
		return "Nonstiff B5";
	}
};

#endif

