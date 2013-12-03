#ifndef NONSTIFF_B1_H
#define NONSTIFF_B1_H

class NonstiffB1 : public BaseIVP {	
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac(0,0) = 2*(1 - y[1]); jac(0,1) = -2*y[0];
		jac(1,0) = y[1];         jac(1,1) = -1 - y[0];
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] =  2*(y[0] - y[0]*y[1]);
		yp[1] =	-1*(y[1] - y[0]*y[1]);
	}
	
public:
	NonstiffB1(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(2);
		_initialCondition[0] = 1;
		_initialCondition[1] = 3;
	}

	virtual const char* GetName() {
		return "Nonstiff B1";
	}
};

#endif

