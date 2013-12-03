#include <core/paramvalue.h>
#include <core/hash.h>

template <class T>
class Allocator {
	typedef T* (*AllocFunction)(Hash<ParamValue>&);
	
	AllocFunction _allocFunction;
	Hash<ParamValue> _params;
	
public:
	Allocator() : _allocFunction(0) {
	}
	
	Allocator(AllocFunction alloc, Hash<ParamValue>& params) {
		_params = params;
		_allocFunction = alloc;
	}
	
	T* Alloc() {
		return _allocFunction(_params);
	}
};
