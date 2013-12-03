#ifndef ALLEN_CAHN_H
#define ALLEN_CAHN_H

class AllenCahn : public TwoSplittingIVP {
	long _n;
	FP _h;
	FP _alpha;
	FP _gamma;

	inline void JacValue(FP* values, long* colInd, long& pos, FP value, long ind) {
		values[pos] = value;
		colInd[pos] = ind;
		pos++;
	}

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		if( split > 2 )
			throw Exception() << "Analytic split " << split << " is not defined.";

		long count = 5*_n*_n - 4*_n;
		FP* values = new FP[count];
		long* colInd = new long[count];
		long* rowPtr = new long[_n*_n];

		FP a = split != 2 ? _alpha/sqr(_h) : 0;
		FP g = split != 1 ? _gamma : 0;

		long pos = 0;
		for( long j = 0; j < _n; j++ ) {
			for( long i = 0; i < _n; i++ ) {
				FP im1jVal = i < _n-1 ? a : 2*a;
				FP ip1jVal = i > 0 ? a : 2*a;
				FP ijm1Val = j < _n-1 ? a : 2*a;
				FP ijp1Val = j > 0 ? a : 2*a;

				rowPtr[j*_n+i] = pos;
				if( j > 0 ) JacValue(values, colInd, pos, ijm1Val, (j-1)*_n+i);
				if( i > 0 ) JacValue(values, colInd, pos, im1jVal, j*_n+i-1);
				JacValue(values, colInd, pos, -4*a + g*(1-3*sqr(y[j*_n+i])), j*_n+i);
				if( i < _n-1 ) JacValue(values, colInd, pos, ip1jVal, j*_n+i+1);
				if( j < _n-1 ) JacValue(values, colInd, pos, ijp1Val, (j+1)*_n+i);
			}
		}

		jac = CSRMat<FP>(values, colInd, rowPtr, _n*_n, _n*_n, count);

		delete [] values;
		delete [] colInd;
		delete [] rowPtr;
	}

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				T yim1j = i > 0 ? y[j*_n+i-1] : y[j*_n+i+1];
				T yip1j = i < _n-1 ? y[j*_n+i+1] : y[j*_n+i-1];
				T yijm1 = j > 0 ? y[(j-1)*_n+i] : y[(j+1)*_n+i];
				T yijp1 = j < _n-1 ? y[(j+1)*_n+i] : y[(j-1)*_n+i];

				yp[j*_n+i] = _alpha*(yim1j + yip1j + yijm1 + yijp1 - 4*y[j*_n+i])/sqr(_h);
			}
		}
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				T yij = y[j*_n+i];
				yp[j*_n+i] = _gamma*yij*(1-sqr(yij));
			}
		}
	}

public:
	AllenCahn(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		SetDefaultFP(params, "tf", 100.);

		_alpha = GetDefaultFP(params, "alpha", 1e-1);
		_gamma = GetDefaultFP(params, "gamma", 1e0);
		_n = GetDefaultLong(params, "N", 20);
		_h = 1./(_n-1);
		
		_initialCondition.Resize(_n*_n);
		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				FP x = i*_h;
				FP y = j*_h;
				_initialCondition[j*_n+i] = 0.4 + 0.1*(x+y) + 0.1*sin(10*x)*sin(20*y);
			}
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("Allen-Cahn")
};

#endif

