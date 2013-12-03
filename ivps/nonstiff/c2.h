#ifndef NONSTIFF_C2_H
#define NONSTIFF_C2_H

class NonstiffC2 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac.Zero();
		jac(0,0) = -1;
		for( long i = 1; i < 9; i++ ) {
			jac(i,i)   = -i;
			jac(i,i-1) =  i-1;
		}
		jac(9,8) = 9;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] =  -y[0];
		yp[1] =   y[0] - 2*y[1];
		yp[2] = 2*y[1] - 3*y[2];
		yp[3] = 3*y[2] - 4*y[3];
		yp[4] = 4*y[3] - 5*y[4];
		yp[5] = 5*y[4] - 6*y[5];
		yp[6] = 6*y[5] - 7*y[6];
		yp[7] = 7*y[6] - 8*y[7];
		yp[8] = 8*y[7] - 9*y[8];
		yp[9] = 9*y[8];
	}

public:
	NonstiffC2(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(10);
		_initialCondition.Zero();
		_initialCondition[0] = 1;
	}
		
	virtual const char* GetName() {
		return "Nonstiff C2";
	}
};

#endif

