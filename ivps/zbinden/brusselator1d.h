#ifndef BRUSSELATOR_1D_H
#define BRUSSELATOR_1D_H

#include <ivps/splitivp.h>

class Brusselator1D : public TwoSplittingIVP {
	long _n;
	FP _alpha;
	FP _coeff;

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		if( split == 0 ) {
			long count = 8*_n - 4;
			FP* values = new FP[count];
			long* colInd = new long[count];
			long* rowPtr = new long[2*_n];

			colInd[0] = 0; values[0] = -2*_coeff + 2*y[0]*y[1] - 4;
			colInd[1] = 1; values[1] = sqr(y[0]);
			colInd[2] = 2; values[2] = _coeff;
			rowPtr[0] = 0;

			colInd[3] = 0; values[3] = 3 - 2*y[0]*y[1];
			colInd[4] = 1; values[4] = -2*_coeff - sqr(y[0]);
			colInd[5] = 3; values[5] = _coeff;
			rowPtr[1] = 3;

			for( long i = 2; i < 2*(_n-1); i += 2 ) {
				FP* val = values + 4*i - 2;
				long* col = colInd + 4*i - 2;

				rowPtr[i] = 6 + 4*(i-2);
				col[0] = i-2; val[0] = _coeff;
				col[1] = i;   val[1] = -2*_coeff + 2*y[i]*y[i+1] - 4;
				col[2] = i+1; val[2] = sqr(y[i]);
				col[3] = i+2; val[3] = _coeff;
				rowPtr[i+1] = 6 + 4*(i-1);
				col[4] = i-1; val[4] = _coeff;
				col[5] = i;   val[5] = 3 - 2*y[i]*y[i+1];
				col[6] = i+1; val[6] = -2*_coeff - sqr(y[i]);
				col[7] = i+3; val[7] = _coeff;
			}

			colInd[count-6] = 2*_n-4;  values[count-6] = _coeff;
			colInd[count-5] = 2*_n-2;  values[count-5] = -2*_coeff + 2*y[2*_n-2]*y[2*_n-1] - 4;
			colInd[count-4] = 2*_n-1;  values[count-4] = sqr(y[2*_n-2]);
			rowPtr[2*_n-2] = 8*_n-10;

			colInd[count-3] = 2*_n-3; values[count-3] = _coeff;
			colInd[count-2] = 2*_n-2; values[count-2] = 3 - 2*y[2*_n-2]*y[2*_n-1];
			colInd[count-1] = 2*_n-1; values[count-1] = -2*_coeff - sqr(y[2*_n-2]);
			rowPtr[2*_n-1] = 8*_n-7;

			jac = CSRMat<FP>(values, colInd, rowPtr, 2*_n, 2*_n, count);

			delete [] values;
			delete [] colInd;
			delete [] rowPtr;
		} else if( split == 1 ) {
			long count = 6*_n - 4;
			FP* values = new FP[count];
			long* colInd = new long[count];
			long* rowPtr = new long[2*_n];

			colInd[0] = 0; values[0] = -2*_coeff;
			colInd[1] = 2; values[1] = _coeff;
			rowPtr[0] = 0;

			colInd[2] = 1; values[2] = -2*_coeff;
			colInd[3] = 3; values[3] = _coeff;
			rowPtr[1] = 2;

			for( long i = 2; i < 2*(_n-1); i += 2 ) {
				FP* val = values + 3*i - 2;
				long* col = colInd + 3*i - 2;

				col[0] = i-2; val[0] = _coeff;
				col[1] = i;   val[1] = -2*_coeff;
				col[2] = i+2; val[2] = _coeff;
				rowPtr[i] = 4 + 3*(i-2);

				col[3] = i-1; val[3] = _coeff;
				col[4] = i+1; val[4] = -2*_coeff;
				col[5] = i+3; val[5] = _coeff;
				rowPtr[i+1] = 4 + 3*(i-1);
			}

			colInd[count-4] = 2*_n-4;  values[count-4] = _coeff;
			colInd[count-3] = 2*_n-2;  values[count-3] = -2*_coeff;
			rowPtr[2*_n-2] = 6*_n-8;

			colInd[count-2] = 2*_n-3; values[count-2] = _coeff;
			colInd[count-1] = 2*_n-1; values[count-1] = -2*_coeff;
			rowPtr[2*_n-1] = 6*_n-6;

			jac = CSRMat<FP>(values, colInd, rowPtr, 2*_n, 2*_n, count);

			delete [] values;
			delete [] colInd;
			delete [] rowPtr;
		} else {
			throw Exception() << "Analytic split " << split << " is not defined.";
		}
	}
	
	template<class T>
	inline void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		long s = y.Size();

		yp[0] = _coeff*(1.-2*y[0]+y[2]);
		yp[1] = _coeff*(3.-2*y[1]+y[3]);

		for( long i = 2; i < s-2; i += 2) {
			yp[i]   = _coeff*(y[i-2] - 2*y[i] + y[i+2]);
			yp[i+1] = _coeff*(y[i-1] - 2*y[i+1] + y[i+3]);
		}

		yp[s-2] = _coeff*(y[s-4] - 2*y[s-2] + 1.);
		yp[s-1] = _coeff*(y[s-3] - 2*y[s-1] + 3.);
	}

	template<class T>	
	inline void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < y.Size(); i += 2 ) {
			yp[i]   = 1. + sqr(y[i])*y[i+1] - 4*y[i];
			yp[i+1] = 3.*y[i] - sqr(y[i])*y[i+1];
		}
	}

public:
	Brusselator1D(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		SetDefaultFP(params, "tf", 10.);
		
		_alpha = GetDefaultFP(params, "alpha", 2e-2);
		_n = GetDefaultLong(params, "N", 500);
		_coeff = _alpha*sqr(_n+1);
		
		_initialCondition.Resize(2*_n);
		for( long i = 0; i < _n; i++ ) {
			_initialCondition[2*i]   = 1 + sin(2*M_PI*(i+1)/(_n+1));
			_initialCondition[2*i+1] = 3;
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("1D Brusselator")
};

#endif

