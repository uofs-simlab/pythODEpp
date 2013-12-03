#ifndef NONSTIFF_D5_H
#define NONSTIFF_D5_H

class NonstiffD5 : public NonstiffD1 {
public:
	NonstiffD5(Hash<ParamValue>& params) : NonstiffD1(params) {
		_initialCondition[0] = 0.1;
		_initialCondition[3] = 4.358898943540673;
	}
		
	virtual const char* GetName() {
		return "Nonstiff D5";
	}
};

#endif

