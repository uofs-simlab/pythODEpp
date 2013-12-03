#ifndef NONSTIFF_C4_H
#define NONSTIFF_C4_H

class NonstiffC4 : public BaseIVP {
protected:
	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
		jac.Zero();
		jac(0,0) = -2;
		jac(0,1) =  1;
		for( long i = 1; i < 50; i++ ) {
			jac(i,i-1) =  1;
			jac(i,i)   = -2;
			jac(i,i+1) =  1;
		}
		jac(50,49) =  1;
		jac(50,50) = -2;
	}   

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = -2*y[0] + y[1];
		for( long i = 1; i < 50; i++ )
			yp[i] = y[i-1] - 2*y[i] + y[i+1];
		yp[50] = y[49] - 2*y[50];	
	}

public:
	NonstiffC4(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(51);
		_initialCondition.Zero();
		_initialCondition[0] = 1;
	}
	
	virtual const char* GetName() {
		return "Nonstiff C4";
	}
};

#endif

