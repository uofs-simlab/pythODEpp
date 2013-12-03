#ifndef NONSTIFF_B4_H
#define NONSTIFF_B4_H

class NonstiffB4 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		FP sumsq = y[0]*y[0] + y[1]*y[1];
		FP root  = sqrt(sumsq);
		FP p23   = pow(sumsq, 2./3);

		jac(0,0) = -y[1]*y[1]*y[2]/p23;
		jac(0,1) = y[0]*y[1]*y[2]/p23 - 1;
		jac(0,2) = -y[0]/root;

		jac(1,0) = y[0]*y[1]*y[2]/p23 + 1;
		jac(1,1) = -y[0]*y[0]*y[2]/p23;
		jac(1,2) = -y[1]/root;

		jac(2,0) = y[1]*y[1]/p23;
		jac(2,1) = -y[0]*y[1]/p23;
		jac(2,2) = 0;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		FP root = sqrt(y[0]*y[0] + y[1]*y[1]);
		yp[0] = -y[1] - y[0]*y[2]/root;
		yp[1] =  y[0] - y[1]*y[2]/root;
		yp[2] =  y[0]/root;
	}

public:
	NonstiffB4(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(3);
		_initialCondition[0] = 3;
		_initialCondition[1] = 0;
		_initialCondition[2] = 0;
	}
	
	virtual const char* GetName() {
		return "Nonstiff B4";
	}
};

#endif

