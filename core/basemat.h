#ifndef BASE_MAT_H
#define BASE_MAT_H

#include <core/common.h>
#include <core/ioobject.h>
#include <core/vec.h>

template <class T>
class BaseMat : public IOObject {
private:
	// Specifically forbid assignment and copy constructor
	BaseMat(const BaseMat<T>& bm) { }
	BaseMat& operator=(const BaseMat& bm) { }

protected:
	long _m, _n;
	T* _elements;

public:
	BaseMat() : _m(0), _n(0), _elements(0) {
#ifdef DEBUGBUILD
		s_Allocations++;
#endif
	}

	virtual ~BaseMat() {
		if( _elements )
			delete [] _elements;
#ifdef DEBUGBUILD
		s_Deletes++;
#endif
	}

	inline long M() const { return _m; }
	inline long N() const { return _n; }

	// Common operations
	virtual void VectorMult(const Vec<T>& vec, Vec<T>& res) const = 0;
	virtual void Factor() = 0;
	virtual void Solve(Vec<T>& b, Vec<T>& x) = 0;

	// Access operators
	const T operator[](long i) const { return _elements[i]; }
	T& operator[](long i) { return _elements[i]; }

#ifdef DEBUGBUILD
protected:		
	static long s_Allocations;
	static long s_Copies;
	static long s_Deletes;

public:
	static long GetAllocations() { return s_Allocations; }
	static long GetCopies() { return s_Copies; }
	static long GetLeaks() { return s_Allocations - s_Deletes; }
	
	static void PrintStats() {
		std::cout << "Matrix allocations: " << GetAllocations() << std::endl;
		std::cout << "Matrix copies:      " << GetCopies() << std::endl;
		std::cout << "Matrix leaks:       " << GetLeaks() << std::endl;
	}
#endif
};

#ifdef DEBUGBUILD
	template<class T> long BaseMat<T>::s_Allocations = 0;
	template<class T> long BaseMat<T>::s_Copies = 0;
	template<class T> long BaseMat<T>::s_Deletes = 0;
#endif

#endif

