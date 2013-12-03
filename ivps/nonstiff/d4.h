#ifndef NONSTIFF_D4_H
#define NONSTIFF_D4_H

class NonstiffD4 : public NonstiffD1 {
public:
	NonstiffD4(Hash<ParamValue>& params) : NonstiffD1(params) {
		_initialCondition[0] = 0.3;
		_initialCondition[3] = 2.3804761428476167;
	}
		
	virtual const char* GetName() {
		return "Nonstiff D4";
	}
};

#endif

