#ifndef NONSTIFF_C1_H
#define NONSTIFF_C1_H

class NonstiffC1 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac.Zero();
		jac(0,0) = -1;
		for( long i = 1; i < 9; i++ ) {
			jac(i,i)   = -1;
			jac(i,i-1) =  1;
		}
		jac(9,8) = 1;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = -y[0];
		yp[1] = y[0] - y[1];
		yp[2] = y[1] - y[2];
		yp[3] = y[2] - y[3];
		yp[4] = y[3] - y[4];
		yp[5] = y[4] - y[5];
		yp[6] = y[5] - y[6];
		yp[7] = y[6] - y[7];
		yp[8] = y[7] - y[8];
		yp[9] = y[8];
	}

public:
	NonstiffC1(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(10);
		_initialCondition.Zero();
		_initialCondition[0] = 1;
	}
		
	virtual const char* GetName() {
		return "Nonstiff C1";
	}
};

#endif

