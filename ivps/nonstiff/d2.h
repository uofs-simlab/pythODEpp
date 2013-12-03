#ifndef NONSTIFF_D2_H
#define NONSTIFF_D2_H

class NonstiffD2 : public NonstiffD1 {
public:
	NonstiffD2(Hash<ParamValue>& params) : NonstiffD1(params) {
		_initialCondition[0] = 0.7;
		_initialCondition[3] = 1.362770287738494;
	}
		
	virtual const char* GetName() {
		return "Nonstiff D2";
	}
};

#endif

