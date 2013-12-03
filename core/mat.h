#ifndef MAT_H
#define MAT_H

#include <core/basemat.h>

template <class T>
class Mat : public BaseMat<T> {
protected:
	Mat<T>* _LU;
	Mat<T>* _P;

public:
	Mat() : _LU(0), _P(0) { }

	Mat(long m, long n) : _LU(0), _P(0) {
		Resize(m,n);
	}

	// Specific definition for the copy constructor	
	Mat(const Mat<T>& m) : _LU(0), _P(0) {
		CopyConstructor(m);
	}   

	// Allow copying from matrices of different types
	template <class U>
	Mat(const Mat<U>& m) : _LU(0), _P(0) {
		CopyConstructor(m);
	}

	virtual ~Mat() {
		if( _LU ) delete _LU;
		if( _P ) delete _P;
	} 

private:
    template <class U>
    inline void CopyConstructor(const Mat<U>& m) {
		this->_m = m.M();
		this->_n = m.N();
		this->_elements = new T[this->_m*this->_n];

		for( long i = 0; i < this->_m; i++ )
			for( long j = 0; j < this->_n; j++ )
				(*this)(i,j) = m(i,j);
#ifdef DEBUGBUILD
		BaseMat<T>::s_Copies++;
#endif
	}

public:
	void Zero() {
		for( long i = 0; i < this->_m*this->_n; i++ )
			this->_elements[i] = 0;
	}

	long NonZeroCount() const {
		long count = 0;
		for( long i = 0; i < this->_m*this->_n; i++ )
			if( this->_elements[i] != 0 ) count++;
		return count;
	}

	void Resize(long m, long n) {
		if( this->_elements ) {
			if( this->_m != m || this->_n != n ) {
				T* old = this->_elements;
				long om = this->_m;
				long on = this->_n;

				this->_m = m; this->_n = n;
				this->_elements = new T[this->_m*this->_n];
				for( long i = 0; i < this->_m && i < om; i++ ) {
					for( long j = 0; j < this->_n && j < on; j++ ) {
						this->_elements[this->_n*i + j] = old[this->_n*i + j];
					}
				}

				delete [] old;
			}
		} else {
			this->_m = m;
			this->_n = n;
			this->_elements = new T[this->_m*this->_n];
		}
	}

	Mat<T>& operator=(const Mat<T>& m) {
		Resize(m._m, m._n);
		for( long i = 0; i < this->_m*this->_n; i++ )
			this->_elements[i] = m._elements[i];
		
		return *this;
	}
	
	// ARITHMETIC OPERATIONS
	Mat<T> operator+(const Mat<T>& m) const {
		return Mat<T>(*this) += m;
	}
	Mat<T>& operator+=(const Mat<T>& m) {
#ifdef DEBUGBUILD
		if( m._m != this->_m && m._n != this->_n )
			std::cerr << "Warning: addition between matrices of incompatible size." << std::endl;
#endif
		for( long i = 0; i < this->_m*this->_n; i++ )
			this->_elements[i] += m._elements[i];
		return *this;
	}

	Mat<T> operator-(const Mat<T>& m) const {
		return Mat<T>(*this) -= m;
	}
	Mat<T>& operator-=(const Mat<T>& m) {
#ifdef DEBUGBUILD
		if( m._m != this->_m && m._n != this->_n )
			std::cerr << "Warning: subtraction between matrices of incompatible size." << std::endl;
#endif
		for( long i = 0; i < this->_m*this->_n; i++ )
			this->_elements[i] -= m._elements[i];
		return *this;
	}

	Mat<T> operator*(const T& v) const {
		return Mat<T>(*this) *= v;
	}
	Mat<T>& operator*=(const T& v) {
		for( long i = 0; i < this->_m*this->_n; i++ ) {
			this->_elements[i] *= v;
		}
		return *this;
	}

	Mat<T> operator/(const T& v) const {
		return Mat<T>(*this) /= v;
	}
	Mat<T>& operator/=(const T& v) {
		for( long i = 0; i < this->_m*this->_n; i++ )
			this->_elements[i] /= v;
		return *this;
	}


	Mat<T> operator-() const {
		Mat<T> m(this->_m,this->_n);
		for( long i = 0; i < this->_m*this->_n; i++ )
			m._elements[i] = -this->_elements[i];
		return m;
	}
			
	Mat<T> Transpose() const {
		Mat<T> m(this->_m,this->_n);
		for( long i = 0; i < this->_m; i++ )
			for( long j = 0; j < this->_n; j++ )
				m(i,j) = (*this)(j,i);
		return m;
	}

	Mat<T> operator*(const Mat<T>& m) const {
#ifdef DEBUGBUILD
		if( m._m != this->_n )
			std::cerr << "Warning: multiplication between matrices of incompatible size." << std::endl;
#endif
		Mat<T> ret(this->_m,this->_n);
		ret.Zero();

		for( int i = 0; i < this->_m; i++ )
			for( int j = 0; j < this->_n; j++ )
				for( int l = 0; l < this->_m; l++ )
					ret(i,j) += (*this)(i,l) * m(l,j);
		return ret;
	}

	virtual void VectorMult(const Vec<T>& vec, Vec<T>& res) const {
#ifdef DEBUGBUILD
		if( vec.Size() != this->_n )
			std::cerr << "Warning: multiplication between matrix and vector of incompatible sizes." << std::endl;
#endif

		for( int i = 0; i < this->_m; i++ ) {
			T sum = 0;
			for( int j = 0; j < this->_n; j++ )
				sum += (*this)(i,j) * vec(j);
			res(i) = sum;
		}
	}
	
	Vec<T> operator*(const Vec<T>& v) const {
		Vec<T> ret(this->_m);
		VectorMult(v,ret);
		return ret;
	}
			
	// Access operators	
	const T operator()(long i, long j) const { return this->_elements[this->_n*i + j]; }
	T& operator()(long i, long j) { return this->_elements[this->_n*i + j]; }

			
	const T MaxNorm() const {
		T mn = 0;
		for( long i = 0; i < this->_m*this->_n; i++ ) {
			T t = fabs(this->_elements[i]);
			if( t > mn ) mn = t;
		}
		return mn;
	}
			
	// Input/ Output
	virtual void Print(std::ostream &out) const {
		for( long i = 0; i < this->_m; i++ ) {
			for( long j = 0; j < this->_n; j++ ) {
				out << (*this)(i,j);
				if( j != this->_n-1 ) out << " ";
			}
			out << std::endl;
		}
	}
	
	virtual void PrintMatlab(std::ostream &out, int printWidth = 0) const {
		out << "[" << std::endl;
		for( long i = 0; i < this->_m; i++ ) {
			for( long j = 0; j < this->_n; j++ ) {
				if( printWidth ) {
					out << " ";
					out.width(printWidth);
					out << (*this)(i,j);
				} else
					out << " " << (*this)(i,j);
			}
			out << std::endl;
		}
		out << "]" << std::endl;
	}

	virtual void Dump(std::ostream &out) const {
		out.write((char*)&this->_m, sizeof(this->_m));
		out.write((char*)&this->_n, sizeof(this->_n));
		out.write((char*)this->_elements, sizeof(T)*this->_m*this->_n);
	}

	virtual void Load(std::istream &in) {
		in.read((char*)&this->_m, sizeof(this->_m));
		in.read((char*)&this->_n, sizeof(this->_n));
		
		if( this->_elements )
			delete [] this->_elements;
		this->_elements = new T[this->_m*this->_n];

		in.read((char*)this->_elements, sizeof(T)*this->_m*this->_n);
	}
		
private:
	inline long SelectPivotRow(long k) const {
		long max = k;
		for( long i = k+1; i < this->_m; i++ )
			if( fabs((*this)(i,k)) > fabs((*this)(max,k)) )
				max = i;
		return max;
	}

	inline void ExchangeRows(long c1, long c2, long r1, long r2) {
		for( long j = c1; j < c2; j++ ) {
			T temp = (*this)(r1,j);
			(*this)(r1,j) = (*this)(r2,j);
			(*this)(r2,j) = temp;
		}
	}

public:
	Mat<T> DestructiveLU() {
		Vec<long> pivot(this->_m);

		long i, j, k;
		T *p_k, *p_row, *p_col;
		FP max;
		T *A = this->_elements;
		
		for (k = 0, p_k = A; k < this->_m; p_k += this->_m, k++) {
			pivot[k] = k;
			max = fabs( *(p_k + k) );
			for (j = k + 1, p_row = p_k + this->_m; j < this->_m; j++, p_row += this->_m) {
				FP temp = fabs(*(p_row + k));
				if ( max < temp) {
					max = temp;
					pivot[k] = j;
					p_col = p_row;
				}
			}

			if (pivot[k] != k)
				for (j = 0; j < this->_m; j++) {
					T swap = *(p_k + j);
					*(p_k + j) = *(p_col + j);
					*(p_col + j) = swap;
				}

			for (i = k+1, p_row = p_k + this->_m; i < this->_m; p_row += this->_m, i++) {
				*(p_row + k) /= *(p_k + k);
			}  

			for (i = k+1, p_row = p_k + this->_m; i < this->_m; p_row += this->_m, i++)
				for (j = k+1; j < this->_m; j++)
					*(p_row + j) -= *(p_row + k) * *(p_k + j);

		}

		Mat<T> pm = Eye(this->_m);
		for( i = 0; i < this->_m; i++ )
			pm.ExchangeRows(0,this->_m,i,pivot[i]);
	
	   return pm;
	}

	virtual void Factor() {
		if( _LU ) delete _LU;
		if( _P ) delete _P;

		_LU = new Mat<T>(*this);
		_P = new Mat<T>(this->_m, this->_n);
		*_P = _LU->DestructiveLU();
	}
	
	void CalcLU(Mat<T>& LU, Mat<T>& P) {
		Factor();
		LU = *_LU;
		P = *_P;
	}

	virtual void Solve(Vec<T>& b, Vec<T>& x) {
		if( !_LU )
			throw Exception() << "Attempted sparse solve without factorizing first\n";

		x = SolveLU(*_LU, *_P, b);
	}
	
	static Vec<T> SolveLU(const Mat<T>& LU, const Mat<T>& P, const Vec<T>& b) {
		Vec<T> x(b.Size());
		Vec<T> pb(P*b);
		
		// Forward substitution
		for( long i = 0; i < b.Size(); i++ ) {
			x(i) = pb(i);
			for( long j = 0; j < i; j++ )
				x(i) -= x(j)*LU(i,j);
		}
		
		// Back substitution
		for( long i = b.Size()-1; i >= 0; i-- ) {
			for( long j = i+1; j < b.Size(); j++ )
				x(i) -= x(j)*LU(i,j);
			x(i) /= LU(i,i);
		}
		
		return x;
	}
			
	static Mat<T> Eye(long s) {
		Mat<T> eye(s,s);
		for( long i = 0; i < s; i++)
			for( long j = 0; j < s; j++)
				eye(i,j) = i == j ? 1 : 0;
		return eye;
	}

	static Mat<T> Rand(long m, long n) {
		Mat<T> r(m,n);
		for( long i = 0; i < m; i++ )
			for( long j = 0; j < n; j++ )
				r(i,j) = 2*(T(rand())/RAND_MAX-(T)0.5);
		return r;
	}
};

template <class T>
Mat<T> operator*(const T& v, const Mat<T>& m) { return m*v; }

#endif

