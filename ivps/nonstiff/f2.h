#ifndef NONSTIFF_F2_H
#define NONSTIFF_F2_H

class NonstiffF2 : public BaseIVP {
protected:
	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp[0] = 55 - (int(t)%2==0?3:1)*y[0]/2;
	}

public:
	NonstiffF2(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);

		_initialCondition.Resize(1);
		_initialCondition[0] = 110;
	}
		
	virtual const char* GetName() {
		return "Nonstiff F2";
	}
};

#endif

