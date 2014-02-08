#ifndef CSR_MAT_H
#define CSR_MAT_H

#include <core/common.h>
#include <core/ioobject.h>
#include <core/vec.h>
#include <core/basemat.h>
#include <core/mat.h>
#include <core/timer.h>

#if USE_SUITESPARSE
	#ifdef __APPLE__
		#include <umfpack.h>
	#else
		#include <suitesparse/umfpack.h>
	#endif
#endif

template <class T>
class CSRMat : public BaseMat<T> {
protected:
	int* _colInd;
	int* _rowPtr;
public:
	void* _factor;
	void* _symbolic;
	long _count;

	void AllocData() {
		this->_elements = new T[_count];
		_colInd = new int[_count];
		_rowPtr = new int[this->_n+1];
	}

public:
	CSRMat() {
		_colInd = 0;
		_rowPtr = 0;
		_factor = 0;
		_symbolic = 0;
	}

	CSRMat(FP* elements, long* col, long* row, long m, long n, long count) {
		_count = count;
		this->_m = m;
		this->_n = n;
		_factor = 0;
		_symbolic = 0;

		AllocData();

		for( long i = 0; i < _count; i++ ) {
			this->_elements[i] = elements[i];
			_colInd[i] = col[i];
		}

		for( long i = 0; i < this->_n; i++ )
			_rowPtr[i] = row[i];
		_rowPtr[this->_n] = _count;
	}

	CSRMat(const Mat<T>& mat) {
		_factor = 0;
		_symbolic = 0;
		FromDense(mat);
	}

	CSRMat(const CSRMat<T>& mat) {
		_factor = 0;
		_symbolic = 0;
		FromSparse(mat);
	}
	
	template <class U>
	CSRMat(const CSRMat<U>& mat) {
		_factor = 0;
		_symbolic = 0;
		FromSparse(mat);
	}

	const long Count() const {
		return _count;
	}

	const int* ColInd() const {
		return _colInd;
	}

	const int* RowPtr() const {
		return _rowPtr;
	}

	virtual ~CSRMat() {
		if( _colInd ) delete [] _colInd;
		if( _rowPtr ) delete [] _rowPtr;
		FreeUMFPack();
	}

	void FreeUMFPack();

	CSRMat<T>& operator=(const Mat<T>& mat) {
		FromDense(mat);
		return *this;
	}

	CSRMat<T>& operator=(const CSRMat<T>& mat) {
		FromSparse(mat);
		return *this;
	}
	
	CSRMat<T>& operator+=(const CSRMat<T>& m) {
#ifdef DEBUGBUILD
		if( m._m != this->_m && m._n != this->_n )
			std::cerr << "Warning: addition between matrices of incompatible size." << std::endl;
#endif
		_count = 0;

		// Loop through once to find how many elements are in the new matrix
		for( long i = 0; i < this->_n; i++ ) {
			long j1 = _rowPtr[i];
			long j2 = m._rowPtr[i];

			for( ; j1 < _rowPtr[i+1] && j2 < m._rowPtr[i+1]; _count++ ) {
				long c1 = _colInd[j1];
				long c2 = m._colInd[j2];
				if( c1 > c2 ) {
					j2++;
				} else if( c2 > c1 ) {
					j1++;
				} else {
					j1++;
					j2++;
				}
			}

			for( ; j1 < _rowPtr[i+1]; _count++, j1++ );
			for( ; j2 < m._rowPtr[i+1]; _count++, j2++ );
		}

		// Allocate and fill everything in now
		T* elementsOld = this->_elements;
		int* rowPtrOld = _rowPtr;
		int* colIndOld = _colInd;
		AllocData();

		long index = 0;
		for( long i = 0; i < this->_n; i++ ) {
			long j1 = rowPtrOld[i];
			long j2 = m._rowPtr[i];
			_rowPtr[i] = index;

			for( ; j1 < rowPtrOld[i+1] && j2 < m._rowPtr[i+1]; index++ ) {
				long c1 = colIndOld[j1];
				long c2 = m._colInd[j2];

				if( c1 > c2 ) {
					_colInd[index]   = c2;
					this->_elements[index] = m._elements[j2];
					j2++;
				} else if( c2 > c1 ) {
					_colInd[index]   = c1;
					this->_elements[index] = elementsOld[j1];
					j1++;
				} else {
					_colInd[index]   = c2;
					this->_elements[index] = elementsOld[j1] + m._elements[j2];
					j1++;
					j2++;
				}
			}

			for( ; j1 < rowPtrOld[i+1]; index++, j1++ ) {
				long c1 = colIndOld[j1];
				_colInd[index]   = c1;
				this->_elements[index] = elementsOld[j1];
			}

			for( ; j2 < m._rowPtr[i+1]; index++, j2++ ) {
				long c2 = m._colInd[j2];
				_colInd[index]   = c2;
				this->_elements[index] = m._elements[j2];
			}
		}
		_rowPtr[this->_n] = _count;

		delete [] elementsOld;
		delete [] rowPtrOld;
		delete [] colIndOld;

		return *this;		
	}

	CSRMat<T> operator+(const CSRMat<T>& m) {
		return CSRMat(*this) += m;
	}
	
	CSRMat<T>& operator-=(const CSRMat<T>& m) {
		return (*this) += m*(-1);
	}
	
	CSRMat<T> operator-(const CSRMat<T>& m) {
		return CSRMat(*this) -= m;
	}

	CSRMat<T> operator*(const T& v) const {
		return CSRMat<T>(*this) *= v;
	}

	CSRMat<T>& operator*=(const T& v) {
		for( long i = 0; i < _count; i++ )
			this->_elements[i] *= v;
		return *this;
	}
	
	CSRMat<T> operator/(const T& v) const {
		return CSRMat<T>(*this) /= v;
	}

	CSRMat<T>& operator/=(const T& v) {
		for( long i = 0; i < _count; i++ )
			this->_elements[i] /= v;
		return *this;
	}

	template <class U>
	void FromSparse(const CSRMat<U>& mat) {
		_count = mat.Count();
		this->_m = mat.M();
		this->_n = mat.N();
		AllocData();

		for( long i = 0; i < _count; i++ ) {
			this->_elements[i] = mat[i];
			_colInd[i] = mat.ColInd()[i];
		}

		for( long i = 0; i <= this->_n; i++ )
			_rowPtr[i] = mat.RowPtr()[i];
	}

	void FromDense(const Mat<T>& mat) {
		_count = mat.NonZeroCount();
		this->_m = mat.M();
		this->_n = mat.N();
		_factor = 0;
		AllocData();

		for( long i = 0, counter = 0; i < this->_n; i++ ) {
			_rowPtr[i] = counter;

			for( long j = 0; j < this->_m; j++ ) {
				FP val = mat(i,j);
				if( val != 0 ) {
					this->_elements[counter] = val;
					_colInd[counter] = j;
					counter++;
				}
			}
		}

		_rowPtr[this->_n] = _count;
	}

	Mat<T> ToDense() {
		Mat<T> mat(this->_m, this->_n);
		mat.Zero();

		for( long i = 0; i < this->_n; i++ )
			for( long j = _rowPtr[i]; j < _rowPtr[i+1]; j++ )
				mat(i,_colInd[j]) = this->_elements[j];
		return mat;
	}

	virtual void VectorMult(const Vec<T>& vec, Vec<T>& res) const {
		for( long i = 0; i < this->_n; i++ ) {
			T sum = 0;
			for( long j = _rowPtr[i]; j < _rowPtr[i+1]; j++ )
				sum += this->_elements[j]*vec[_colInd[j]];
			res[i] = sum;
		}
	}

	virtual void Factor();
	virtual void Solve(Vec<T>& b, Vec<T>& x);

	virtual void Dump(std::ostream &out) const {
	}

	virtual void Load(std::istream &in) {
	}

	virtual void Print(std::ostream &out, int printWidth = 0) const {
		out << "\n";
		for( long i = 0; i < _count; i++ ) {
			if( printWidth ) {
				out << " ";
				out.width(printWidth);
				out  << this->_elements[i];
			} else
				out << " " << this->_elements[i];
		}
		out << "\n";

		for( long i = 0; i < _count; i++ ) {
			if( printWidth ) {
				out << " ";
				out.width(printWidth);
				out  << _colInd[i];
			} else
				out << " " << _colInd[i];
		}
		out << "\n";

		for( long i = 0; i < this->_n; i++ ) {
			if( printWidth ) {
				out << " ";
				out.width(printWidth);
				out  << _rowPtr[i];
			} else
				out << " " << _rowPtr[i];
		}
		out << "\n";
	}

	virtual void PrintMatlab(std::ostream &out, int printWidth = 0 ) const {
		out << "m = " << this->_m << ";\n";
		out << "n = " << this->_n << ";\n";
		out << "i = [";
		for( long i = 0; i < this->_n; i++ ) {
			for( long j = 0; j < _rowPtr[i+1]-_rowPtr[i]; j++ )
				out << " " << i+1;
		}
		out << "];\n";
		
		out << "j = [";
		for( long i = 0; i < _count; i++ )
			out << " " << _colInd[i]+1;
		out << "];\n";

		out << "s = [";
		for( long i = 0; i < _count; i++ ) {
			out << " ";
			if( printWidth )
				out.width(printWidth);
			out << this->_elements[i];
		}
		out << "];\n";
	}
	
	virtual void PrintPython(std::ostream &out, int printWidth = 0 ) const {
		out << "m = " << this->_m << ";\n";
		out << "n = " << this->_n << ";\n";
		out << "i = [";
		for( long i = 0; i < this->_n; i++ )
			out << (i == 0 ? " " : ", ") << _rowPtr[i];
		out << ", " << _rowPtr[this->_n] << "];\n";
		
		out << "j = [";
		for( long i = 0; i < _count; i++ )
			out << (i == 0 ? "" : ",") << " " << _colInd[i];
		out << "];\n";

		out << "s = [";
		for( long i = 0; i < _count; i++ ) {
			out << (i == 0 ? "" : ",") << " ";
			if( printWidth )
				out.width(printWidth);
			out << this->_elements[i];
		}
		out << "];\n";
	}

	static CSRMat<T> Eye(long s) {
		long* rows = new long[s];
		long* cols = new long[s];
		FP* data = new FP[s];

		for( long i = 0; i < s; i++ ) {
			rows[i] = i;
			cols[i] = i;
			data[i] = 1;
		}

		CSRMat eyeMat(data, cols, rows, s, s, s);

		delete [] data;
		delete [] cols;
		delete [] rows;

		return eyeMat;
	}
};

template<> inline
void CSRMat<FP>::FreeUMFPack() {
#ifdef USE_SUITESPARSE
	if( _factor )
		umfpack_di_free_numeric(&_factor);
	if( _symbolic )
		umfpack_di_free_symbolic(&_symbolic);
#else
	throw Exception() << "Not compiled with support for sparse solvers.";
#endif
}

template<> inline
void CSRMat<CFP>::FreeUMFPack() {
#ifdef USE_SUITESPARSE
	if( _factor )
		umfpack_zi_free_numeric(&_factor);
	if( _symbolic )
		umfpack_zi_free_symbolic(&_symbolic);
#else
	throw Exception() << "Not compiled with support for sparse solvers.";
#endif
}

// TODO: FIXME: Think of a better way to leave error messages here
template<> inline
void CSRMat<FP>::Factor() {
#ifdef USE_SUITESPARSE
	if( _symbolic )
		umfpack_di_free_symbolic(&_symbolic);
	if( _factor )
		umfpack_di_free_numeric(&_factor);

	double info[UMFPACK_INFO];
	int status;

	//Timer st;
	if( (status = umfpack_di_symbolic(this->_n, this->_m, _rowPtr, _colInd, this->_elements, &_symbolic, 0, info)) < 0 ) {
		umfpack_di_report_info (0, info);
		umfpack_di_report_status (0, status);
		throw Exception() << "umfpack_di_symbolic failed.";
	}
	//printf("    symbolic factor time: %dms\n", (int)st.msec());

	//Timer nf;
	if( (status = umfpack_di_numeric(_rowPtr, _colInd, this->_elements, _symbolic, &_factor, 0, info)) < 0 ) {
		umfpack_di_report_info (0, info);
		umfpack_di_report_status (0, status);
		throw Exception() << "umfpack_di_numeric failed.";
	}
	//printf("    numeric factor time: %dms\n", (int)nf.msec());
#else
	throw Exception() << "Not compiled with support for sparse solvers.";
#endif
}

template<> inline
void CSRMat<CFP>::Factor() {
#ifdef USE_SUITESPARSE
	if( _symbolic )
		umfpack_zi_free_symbolic(&_symbolic);
	if( _factor )
		umfpack_zi_free_numeric(&_factor);

	FP* re = new FP[_count];
	FP* im = new FP[_count];
	for( long i = 0; i < _count; i++ ) {
		re[i] = this->_elements[i].real();
		im[i] = -this->_elements[i].imag();
	}

	umfpack_zi_symbolic(this->_n, this->_m, _rowPtr, _colInd, re, im, &_symbolic, 0, 0);
	umfpack_zi_numeric(_rowPtr, _colInd, re, im, _symbolic, &_factor, 0, 0);

	delete [] re;
	delete [] im;
#else
	throw Exception() << "Not compiled with support for sparse solvers.";
#endif
}

template<> inline 
void CSRMat<FP>::Solve(Vec<FP>& b, Vec<FP>& x) {
#ifdef USE_SUITESPARSE
	if( !_factor )
		throw Exception() << "Attempted sparse solve without factorizing first\n";
	umfpack_di_solve(UMFPACK_At, _rowPtr, _colInd, this->_elements, *x, *b, _factor, 0, 0);
#else
	throw Exception() << "Not compiled with support for sparse solvers.";
#endif
}

template<> inline 
void CSRMat<CFP>::Solve(Vec<CFP>& b, Vec<CFP>& x) {
#ifdef USE_SUITESPARSE
	if( !_factor )
		throw Exception() << "Attempted sparse solve without factorizing first\n";
	Vec<FP> rex(x.Size());
	Vec<FP> imx(x.Size());
	Vec<FP> reb = VecReal(b);
	Vec<FP> imb = VecImag(b);
	
	FP* re = new FP[_count];
	FP* im = new FP[_count];
	for( long i = 0; i < _count; i++ ) {
		re[i] = this->_elements[i].real();
		im[i] = -this->_elements[i].imag();
	}

	umfpack_zi_solve(UMFPACK_At, _rowPtr, _colInd, re, im, *rex, *imx, *reb, *imb, _factor, 0, 0);

	for( long i = 0; i < x.Size(); i++ )
		x[i] = CFP(rex[i], imx[i]);
	
	delete [] re;
	delete [] im;
#else
	throw Exception() << "Not compiled with support for sparse solvers.";
#endif
}

template <class T>
CSRMat<T> operator*(const T& v, const CSRMat<T>& m) { return m*v; }

#endif

