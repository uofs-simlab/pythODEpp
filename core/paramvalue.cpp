#include <core/paramvalue.h>

void ParamValue::Clean() {
	if( _strValue ) {
		delete [] _strValue;
	}
	
	_strValue = 0;
	_floatValue = 0;
	_longValue = 0;
}

ParamValue::ParamValue() : _strValue(0), _floatValue(0), _longValue(0) {
}

ParamValue::ParamValue(const char* value) : _strValue(0), _floatValue(0), _longValue(0) {
	SetString(value);
}

ParamValue::ParamValue(FP value) : _strValue(0), _floatValue(0), _longValue(0) {
	SetFP(value);
}

ParamValue::ParamValue(long value) : _strValue(0), _floatValue(0), _longValue(0) {
	SetLong(value);
}

ParamValue::ParamValue(const ParamValue& pv) : _strValue(0),
	_floatValue(pv._floatValue), _longValue(pv._longValue) {
	if( pv._strValue ) {
		long len = strlen(pv._strValue);
		_strValue = new char[len+1];
		strcpy(_strValue, pv._strValue);
    }
}

ParamValue::~ParamValue() {
	if( _strValue ) delete [] _strValue;
}

void ParamValue::SetLong(const long value) {
	if( _floatValue == value )
		return;
	
	Clean();
	_floatValue = value;
	_longValue = long(value);
}

void ParamValue::SetFP(const FP value) {
	if( _floatValue == value )
		return;
	
	Clean();
	_floatValue = value;
	_longValue = long(value);
}

void ParamValue::SetString(const char* str) {
	if( _strValue && strcmp(str, _strValue) == 0 )
		return;
	
	Clean();
	if( !str ) return;
	
	long len = strlen(str);
	_strValue = new char[len+1];
	strcpy(_strValue,str);
	_floatValue = atof(_strValue);
	_longValue = atoi(_strValue);
}

long ParamValue::GetLong() {
	return _longValue;
}

FP ParamValue::GetFP() {
	return _floatValue;
}

const char* ParamValue::GetString() {
	if( !_strValue ) {
		char buf[128];
		long len = snprintf(buf,sizeof(buf),"%f",_floatValue);
		_strValue = new char[len+1];
		strcpy(_strValue,buf);
	}
	return _strValue;
}

ParamValue& ParamValue::operator=(const ParamValue& pv) {
	_strValue   = 0;
	_floatValue = pv._floatValue;
	_longValue  = pv._longValue;

	if( pv._strValue ) {
		long len = strlen(pv._strValue);
		_strValue = new char[len+1];
		strcpy(_strValue, pv._strValue);
    }

	return *this;
}

bool ParamValue::operator<(ParamValue& pv) {
	return strcmp(GetString(), pv.GetString()) < 0;
}

bool ParamValue::operator>(ParamValue& pv) {
	return strcmp(GetString(), pv.GetString()) > 0;
}

bool ParamValue::operator<=(ParamValue& pv) {
	return strcmp(GetString(), pv.GetString()) <= 0;
}

bool ParamValue::operator>=(ParamValue& pv) {
	return strcmp(GetString(), pv.GetString()) >= 0;
}

bool ParamValue::operator==(ParamValue& pv) {
	return strcmp(GetString(), pv.GetString()) == 0;
}

bool ParamValue::operator!=(ParamValue& pv) {
	return strcmp(GetString(), pv.GetString()) != 0;
}

long SetDefaultLong(Hash<ParamValue>& params, const char* name, long d) {
	ParamValue* pv = params.Get(name);
	if( pv ) return pv->GetLong();
	params[name].SetLong(d);
	return d;
}

long GetDefaultLong(Hash<ParamValue>& params, const char* name, long d) {
	ParamValue* pv = params.Get(name);
	if( pv ) return pv->GetLong();
	return d;
}

FP SetDefaultFP(Hash<ParamValue>& params, const char* name, FP d) {
	ParamValue* pv = params.Get(name);
	if( pv ) return pv->GetFP();
	params[name].SetFP(d);
	return d;
}

FP GetDefaultFP(Hash<ParamValue>& params, const char* name, FP d) {
	ParamValue* pv = params.Get(name);
	if( pv ) return pv->GetFP();
	return d;
}

const char* SetDefaultString(Hash<ParamValue>& params, const char* name, const char* d) {
	ParamValue* pv = params.Get(name);
	if( pv ) return pv->GetString();
	params[name].SetString(d);
	return params[name].GetString();
}

const char* GetDefaultString(Hash<ParamValue>& params, const char* name, const char* d) {
	ParamValue* pv = params.Get(name);
	if( pv ) return pv->GetString();
	return params[name].GetString();
}

