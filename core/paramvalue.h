#ifndef PARAM_VALUE_H
#define PARAM_VALUE_H

#include <core/common.h>
#include <core/hash.h>

class ParamValue {
	char* _strValue;
	FP    _floatValue;
	long  _longValue;

	void Clean();
	
public:
	ParamValue();
	ParamValue(const char* value);
	ParamValue(FP value);
	ParamValue(long value);
	ParamValue(const ParamValue& pv);
	~ParamValue();

	void SetLong(const long value);
	void SetFP(const FP value);
	void SetString(const char* rhs);
	
	long GetLong();
	FP GetFP();
	const char* GetString();

	ParamValue& operator=(const ParamValue& pv);
	bool operator<(ParamValue& pv);
	bool operator>(ParamValue& pv);
	bool operator<=(ParamValue& pv);
	bool operator>=(ParamValue& pv);
	bool operator==(ParamValue& pv);
	bool operator!=(ParamValue& pv);
};

long SetDefaultLong(Hash<ParamValue>& params, const char* name, long d);
long GetDefaultLong(Hash<ParamValue>& params, const char* name, long d);
FP SetDefaultFP(Hash<ParamValue>& params, const char* name, FP d);
FP GetDefaultFP(Hash<ParamValue>& params, const char* name, FP d);
const char* SetDefaultString(Hash<ParamValue>& params, const char* name, const char* d);
const char* GetDefaultString(Hash<ParamValue>& params, const char* name, const char* d);

#endif

