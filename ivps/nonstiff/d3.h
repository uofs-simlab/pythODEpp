#ifndef NONSTIFF_D3_H
#define NONSTIFF_D3_H

class NonstiffD3 : public NonstiffD1 {
public:
	NonstiffD3(Hash<ParamValue>& params) : NonstiffD1(params) {
		_initialCondition[0] = 0.5;
		_initialCondition[3] = 1.7320508075688772;
	}
		
	virtual const char* GetName() {
		return "Nonstiff D3";
	}
};

#endif

