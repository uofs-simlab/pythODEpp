#ifndef VDPOL_H
#define VDPOL_H

class VDPOL : public BaseIVP {
protected:
	FP _param;

	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = 0;
		jac(0,1) = 1;
		jac(1,0) = (-2*y(0)*y(1)-1)/_param;
		jac(1,1) = (1-sqr(y(0)))/_param;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp(0) = y(1);
		yp(1) = ((1-sqr(y(0)))*y(1)-y(0))/_param;
	}

public:
	VDPOL(Hash<ParamValue>& params) : BaseIVP(params), _param(1e-6) {
		params["tf"].SetFP(2);

		_initialCondition.Resize(2);
		_initialCondition(0) = 2;
		_initialCondition(1) = -0.66;
	}

	virtual const char* GetName() {
		return "Van der Pol problem";
	}
};

#endif

