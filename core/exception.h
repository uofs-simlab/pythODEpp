#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <core/common.h>

class Exception {
	std::stringstream _message;

public:
	Exception() { }
	Exception(const Exception& e) {
		_message << e._message.rdbuf();
	}
	
	template <class T>
	Exception& operator<<(T v) {
		_message << v;
		return *this;
	}
	 
	operator std::string () {
		return _message.str();
	}

	operator const char* () {
		return _message.str().c_str();
	}
};

#endif

