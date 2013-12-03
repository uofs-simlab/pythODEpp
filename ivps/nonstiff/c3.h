#ifndef NONSTIFF_C3_H
#define NONSTIFF_C3_H

class NonstiffC3 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac.Zero();
		jac(0,0) = -2;
		jac(0,1) =  1;
		for( long i = 1; i < 9; i++ ) {
			jac(i,i-1) =  1;
			jac(i,i)   = -2;
			jac(i,i+1) =  1;
		}
		jac(9,8) =  1;
		jac(9,9) = -2;
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = -2*y[0] + y[1];
		for( long i = 1; i < 9; i++ )
			yp[i] = y[i-1] - 2*y[i] + y[i+1];
		yp[9] = y[8] - 2*y[9];
	}

public:
	NonstiffC3(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(10);
		_initialCondition.Zero();
		_initialCondition[0] = 1;
	}
		
	virtual const char* GetName() {
		return "Nonstiff C3";
	}
};

#endif

