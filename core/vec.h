#pragma once

#include <core/common.h>
#include <core/exception.h>
#include <core/ioobject.h>

#ifdef USE_ADOL_C
	#include <adolc/adolc.h>
	#include <adolc/adolc_sparse.h>
#endif

template <class T>
class Vec : public IOObject {
protected:
	long _size;
	T* _elements;

public:
	Vec() : _size(0), _elements(0) {
#ifdef DEBUGBUILD
		s_Allocations++;
#endif
	}
	
	Vec(long size) : _size(0), _elements(0) {
		Resize(size);
#ifdef DEBUGBUILD
		s_Allocations++;
#endif
	}

private:
	template <class U>
	inline void CopyConstructor(const Vec<U>& sv) {
		_size = sv.Size();
		_elements = new T[_size];
		for( long i = 0; i < _size; i++ )
			_elements[i] = sv[i];
#ifdef DEBUGBUILD
		s_Allocations++;
		s_Copies++;
#endif
	}

public:

	Vec(const Vec<T>& sv) {
		CopyConstructor(sv);
	}

	// Allow copying from vectors of different types
	template <class U>
	Vec(const Vec<U>& sv) {
		CopyConstructor(sv);
	}
	
	virtual ~Vec() {
		if( _elements )
			delete [] _elements;
#ifdef DEBUGBUILD
		s_Deletes++;
#endif
	}
	
	void Resize(long size) {
		if( _elements ) {
			if( _size != size ) {
				T* old = _elements;
				long oldSize = _size;
				
				_size = size;
				_elements = new T[size];
				for( long i = 0; i < _size && i < oldSize; i++ )
					_elements[i] = old[i];
				
				delete [] old;
			}
		} else {
			_size = size;
			_elements = new T[size];
		}
	}
	
	long Size() const {
		return _size;
	}

	void Zero() {
		for( long i = 0; i < _size; i++ )
			_elements[i] = 0;
	}
	
	T Sum() const {
		T sum = 0;
		for( long i = 0; i < _size; i++ )
			sum += _elements[i];
		return sum;
	}

	T Norm() const {
		T norm = 0;
		for( long i = 0; i < _size; i++ )
			norm += _elements[i]*_elements[i];
		return sqrt(norm);
	}

	T InfNorm() const {
		T norm = 0;
		for( long i = 0; i < _size; i++ ) {
			T temp = fabs(_elements[i]);
			if( temp > norm )
				norm = temp;
		}
		return norm;
	}

	T RMS() const {
		T rms = 0;
		for( long i = 0; i < _size; i++ )
			rms += _elements[i]*_elements[i];
		rms /= _size;
		return sqrt(rms);
	}

	Vec<T>& operator=(const Vec<T>& sv) {
		Resize(sv.Size());
		for( long i = 0; i < _size; i++ )
			_elements[i] = sv[i];
		
		return *this;
	}
	
	Vec<T> operator+(const Vec<T>& sv) const { return Vec<T>(*this) += sv; }
	template <class C>
	Vec<T> operator+(const C& c) const { return Vec<T>(*this) += c; }
	Vec<T> operator-(const Vec<T>& sv) const { return Vec<T>(*this) -= sv; }
	template <class C>
	Vec<T> operator-(const C& c) const { return Vec<T>(*this) -= c; }
	Vec<T> operator*(const Vec<T>& sv) const { return Vec<T>(*this) *= sv; }
	Vec<T> operator/(const Vec<T>& sv) const { return Vec<T>(*this) /= sv; }
	template <class C>
	Vec<T> operator*(const C& c) const { return Vec<T>(*this) *= c; }
	template <class C>
	Vec<T> operator/(const C& c) const { return Vec<T>(*this) /= c; }
	Vec<T> operator-() const {
		Vec<T> ret(_size);
		for( long i = 0; i < _size; i++ )
			ret(i) = -(*this)(i);
		return ret;
	}

	Vec<T>& operator+=(const Vec<T>& sv) {
#ifdef DEBUGBUILD
		if( _size != sv.Size() )
			throw Exception() << "Addition between vectors of incompatible sizes.";
#endif
		for( long i = 0; i < _size && i < sv.Size(); i++ )
			_elements[i] += sv[i];
		return *this;
	}
	
	template <class C>
	Vec<T>& operator+=(const C& c) {
		for( long i = 0; i < _size && i < _size; i++ )
			_elements[i] += c;
		return *this;
	}
	
	Vec<T>& operator-=(const Vec<T>& sv) {
#ifdef DEBUGBUILD
		if( _size != sv.Size() )
			throw Exception() << "Subtraction between vectors of incompatible sizes.";
#endif
		for( long i = 0; i < _size && i < sv.Size(); i++ )
			_elements[i] -= sv[i];
		return *this;
	}
	
	template <class C>
	Vec<T>& operator-=(const C& c) {
		for( long i = 0; i < _size && i < _size; i++ )
			_elements[i] -= c;
		return *this;
	}
	
	Vec<T>& operator*=(const Vec<T>& sv) {
#ifdef DEBUGBUILD
		if( _size != sv.Size() )
			throw Exception() << "Multiplication between vectors of incompatible sizes.";
#endif
		for( long i = 0; i < _size && i < sv.Size(); i++ )
			_elements[i] *= sv[i];
		return *this;
	}
	
	Vec<T>& operator/=(const Vec<T>& sv) {
#ifdef DEBUGBUILD
		if( _size != sv.Size() )
			throw Exception() << "Division between vectors of incompatible sizes.";
#endif
		for( long i = 0; i < _size; i++ )
			_elements[i] /= sv[i];
		return *this;
	}
	
	template <class C>
	Vec<T>& operator*=(const C& c) {
		for( long i = 0; i < _size; i++ )
			_elements[i] *= c;
		return *this;
	}
	
	template <class C>
	Vec<T>& operator/=(const C& c) {
		for( long i = 0; i < _size; i++ )
			_elements[i] /= c;
		return *this;
	}

	inline const T operator[](long i) const {
#ifdef DEBUGBUILD
		if( i < 0 || i >= _size )
			throw Exception() << "Index out of bounds.";
#endif
		return _elements[i];
	}

	inline T& operator[](long i) {
#ifdef DEBUGBUILD
		if( i < 0 || i >= _size )
			throw Exception() << "Index out of bounds.";
#endif
		return _elements[i];
	}
	
	inline const T operator()(long i) const {
#ifdef DEBUGBUILD
		if( i < 0 || i >= _size )
			throw Exception() << "Index out of bounds.";
#endif
		return _elements[i];
	}

	inline T& operator()(long i) {
#ifdef DEBUGBUILD
		if( i < 0 || i >= _size )
			throw Exception() << "Index out of bounds.";
#endif
		return _elements[i];
	}

	T* operator*() {
		return _elements;
	}
	
	bool IsNan() const {
		for( long i = 0; i < _size; i++ )
			if( _elements[i] != _elements[i] )
				return true;
		return false;
	}

	Vec<T> Splice(long i1, long i2) const {
		Vec<T> splice(i2-i1);
		for( long i = 0; i < splice.Size(); i++ )
			splice(i) = _elements[i+i1];
		return splice;
	}

	void AssignSplice(const Vec<T>& v, long index) {
		for( long i = 0; i < v.Size(); i++ )
			(*this)(i+index) = v(i);
	}

	Vec<T> Diff() const {
		Vec<T> diff(_size-1);
		for( long i = 0; i < _size-1; i++ )
			diff(i) = _elements[i+1] - _elements[i];
		return diff;
	}

	Vec<T> Power(T e) const {
		Vec<T> ret(_size);
		for( long i = 0; i < _size; i++ )
			ret[i] = pow(_elements[i],e);
		return ret;
	}

	Vec<T> RealPower(T e) const {
		Vec<T> ret(_size);
		for( long i = 0; i < _size; i++ )
			ret[i] = realpow(_elements[i],e);
		return ret;
	}


	Vec<T> Exp() const {
		Vec<T> ret(_size);
		for( long i = 0; i < _size; i++ )
			ret[i] = exp(_elements[i]);
		return ret;
	}
	
	Vec<T> Maximum(T m) const {
		Vec<T> ret(_size);
		for( long i = 0; i < _size; i++ )
			ret[i] = fmax(_elements[i],m);
		return ret;
	}
	
	Vec<T> Minimum(T m) const {
		Vec<T> ret(_size);
		for( long i = 0; i < _size; i++ )
			ret[i] = fmin(_elements[i],m);
		return ret;
	}

	T PeriodicCentralDifference(long index, long deriv, long order, FP dx, long stride = 1) const {
		if( deriv == 1 ) {
			switch( order ) {
			case 2: {
				long indices[] = { (_size+index-stride)%_size,
								   (_size+index+stride)%_size };
				FP coeffs[] = { -1./2, 1./2 };
				T result = 0;
				for( long i = 0; i < 2; i++ )
					result += coeffs[i]*(*this)(indices[i]);
				return result/dx;
			}
			case 4: {
				long indices[] = { (_size+index-2*stride)%_size,
								   (_size+index-1*stride)%_size,
								   (_size+index+1*stride)%_size,
								   (_size+index+2*stride)%_size };
				FP coeffs[] = { 1./12, -2./3, 2./3, -1./12 };
				T result = 0;
				for( long i = 0; i < 4; i++ )
					result += coeffs[i]*(*this)(indices[i]);
				return result/dx;
			}
			case 6: {
				long indices[] = { (_size+index-3*stride)%_size,
								   (_size+index-2*stride)%_size,
								   (_size+index-1*stride)%_size,
								   (_size+index+1*stride)%_size,
								   (_size+index+2*stride)%_size,
								   (_size+index+3*stride)%_size };
				FP coeffs[] = { -1./60, 3./20, -3./4, 3./4, -3./20, 1./60 };
				T result = 0;
				for( long i = 0; i < 6; i++ )
					result += coeffs[i]*(*this)(indices[i]);
				return result/dx;
			} }
		} else if( deriv == 2 ) {
			switch( order ) {
			case 2: {
				long indices[] = { (_size+index-stride)%_size,
								   index,
								   (_size+index+stride)%_size };
				FP coeffs[] = { 1., -2., 1. };
				T result = 0;
				for( long i = 0; i < 3; i++ )
					result += coeffs[i]*(*this)(indices[i]);
				return result/(dx*dx);
			}
			case 4: {
				long indices[] = { (_size+index-2*stride)%_size,
								   (_size+index-1*stride)%_size,
								   index,
								   (_size+index+1*stride)%_size,
								   (_size+index+2*stride)%_size };
				FP coeffs[] = { -1./12, 4./3, -5./2, 4./3, -1./12 };
				T result = 0;
				for( long i = 0; i < 5; i++ )
					result += coeffs[i]*(*this)(indices[i]);
				return result/(dx*dx);
			}
			case 6: {
				long indices[] = { (_size+index-3*stride)%_size,
								   (_size+index-2*stride)%_size,
								   (_size+index-1*stride)%_size,
								   index,
								   (_size+index+1*stride)%_size,
								   (_size+index+2*stride)%_size,
								   (_size+index+3*stride)%_size };
				FP coeffs[] = { 1./90, -3./20, 3./2, -49./18, 3./2, -3./20, 1./90 };
				T result = 0;
				for( long i = 0; i < 7; i++ )
					result += coeffs[i]*(*this)(indices[i]);
				return result/(dx*dx);
			} }
		}

		throw Exception() << "Finite difference derivative " << deriv << ", order " << order << " is unsupported.";
	}

	void AddScaled(T s, const Vec<T>& v) {
#ifdef DEBUGBUILD
		if( _size != v.Size() )
			throw Exception() << "Adding scaled vectors to incompatible size.";
#endif
		for( long i = 0; i < _size; i++ )
			_elements[i] += s*v._elements[i];
	}

	// Input/ Output
	virtual void PrintMatlab(std::ostream &out = std::cout) const {
		out << "[" << std::endl;
		for( long j = 0; j < _size; j++ )
			out << " " << (*this)(j) << std::endl;
		out << "]" << std::endl;
	}
	
	virtual void Print(std::ostream &out = std::cout) const {
		for( long j = 0; j < _size; j++ )
			out << (*this)(j) << " ";
		out << std::endl;
	}

	virtual void Dump(std::ostream &out) const {
		out.write((char*)&_size, sizeof(_size));
		out.write((char*)_elements, sizeof(T)*_size);
	}

	virtual void Load(std::istream &in) {
		in.read((char*)&_size, sizeof(_size));
		if( _elements )
			delete [] _elements;
		_elements = new T[_size];

		in.read((char*)_elements, sizeof(T)*_size);
	}

	static Vec<T> Rand(long s) {
		Vec<T> r(s);
		for( long i = 0; i < s; i++ )
			r(i) = 2*(T(rand())/RAND_MAX-(T)0.5);
		return r;
	}
	
	static Vec<T> Zeros(long s) {
		Vec<T> r(s);
		for( long i = 0; i < s; i++ )
			r(i) = 0;
		return r;
	}

	static Vec<T> Ones(long s) {
		Vec<T> r(s);
		for( long i = 0; i < s; i++ )
			r(i) = 1;
		return r;
	}

#ifdef DEBUGBUILD
private:
	static long s_Allocations;
	static long s_Copies;
	static long s_Deletes;

public:
	static long GetAllocations() { return s_Allocations; }
	static long GetCopies() { return s_Copies; }
	static long GetLeaks() { return s_Allocations - s_Deletes; }
	
	static void PrintStats() {
		std::cout << "Vector allocations: " << GetAllocations() << std::endl;
		std::cout << "Vector copies:      " << GetCopies() << std::endl;
		std::cout << "Vector leaks:       " << GetLeaks() << std::endl;
	}
#endif
};

#ifdef DEBUGBUILD
	template<class T> long Vec<T>::s_Allocations = 0;
	template<class T> long Vec<T>::s_Copies = 0;
	template<class T> long Vec<T>::s_Deletes = 0;
#endif

template <class T, class C>
Vec<T> operator+(const C& c, const Vec<T>& sv) { return sv + c; }

template <class T, class C>
Vec<T> operator-(const C& c, const Vec<T>& sv) { return -(sv - c); }

template <class T, class C>
Vec<T> operator*(const C& c, const Vec<T>& sv) { return sv*c; }

Vec<FP> VecReal(const Vec<CFP>& cv);
Vec<FP> VecImag(const Vec<CFP>& cv);

