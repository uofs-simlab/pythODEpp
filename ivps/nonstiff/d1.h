#ifndef NONSTIFF_D1_H
#define NONSTIFF_D1_H

class NonstiffD1 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		FP y2s = sqr(y[2]);
		FP y3s = sqr(y[3]);
		FP denom = pow(y2s + y3s, -2.5);

		jac.Zero();
		jac(0,2) = 1;
		jac(1,3) = 1;
		jac(2,0) = (2*y2s-y3s)   * denom;
		jac(2,1) = (3*y[2]*y[3]) * denom;
		jac(3,0) = (3*y[2]*y[3]) * denom;
		jac(3,1) = (2*y3s-y2s)   * denom;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		FP denom = pow(sqr(y[0]) + sqr(y[1]), -1.5);

		yp[0] = y[2];
		yp[1] = y[3];
		yp[2] = -y[0]*denom;
		yp[3] = -y[1]*denom;
	}

public:
	NonstiffD1(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(4);
		_initialCondition[0] = 0.9;
		_initialCondition[1] = 0;
		_initialCondition[2] = 0;
		_initialCondition[3] = 1.1055415967851334;
	}
		
	virtual const char* GetName() {
		return "Nonstiff D1";
	}
};

#endif

