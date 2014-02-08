#ifndef ADVECTION_DIFFUSION_1D_H
#define ADVECTION_DIFFUSION_1D_H

#include <ivps/splitivp.h>

class AdvectionDiffusion1D : public TwoSplittingIVP {
	long _n;
	FP _advection;
	FP _diffusivity;

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		FP diff = _diffusivity*_n*_n;
		FP adv  = _advection*_n/2;
 
		if( split == 0 ) {
			long count = 3*_n;
			FP* values = new FP[count];
			long* colInd = new long[count];
			long* rowPtr = new long[_n];

			rowPtr[0] = 0;
			colInd[0] = 0;    values[0] = -2*diff;
			colInd[1] = 1;    values[1] = diff + adv;
			colInd[2] = _n-1; values[2] = diff - adv;

			for( long i = 1; i < _n-1; i++ ) {
				rowPtr[i] = 3*i;
				colInd[3*i]   = i-1; values[3*i]   = diff - adv;
				colInd[3*i+1] = i;   values[3*i+1] = -2*diff;
				colInd[3*i+2] = i+1; values[3*i+2] = diff + adv;
			}

			rowPtr[_n-1] = count-3;
			colInd[count-3] = 0;    values[count-3] = diff + adv;
			colInd[count-2] = _n-2; values[count-2] = diff - adv;
			colInd[count-1] = _n-1; values[count-1] = -2*diff;

			jac = CSRMat<FP>(values, colInd, rowPtr, _n, _n, count);

			delete [] values;
			delete [] colInd;
			delete [] rowPtr;
		} else if( split == 1 ) {
			long count = 3*_n;
			FP* values = new FP[count];
			long* colInd = new long[count];
			long* rowPtr = new long[_n];

			rowPtr[0] = 0;
			colInd[0] = 0;    values[0] = -2*diff;
			colInd[1] = 1;    values[1] = diff;
			colInd[2] = _n-1; values[2] = diff;

			for( long i = 1; i < _n-1; i++ ) {
				rowPtr[i] = 3*i;
				colInd[3*i]   = i-1; values[3*i]   = diff;
				colInd[3*i+1] = i;   values[3*i+1] = -2*diff;
				colInd[3*i+2] = i+1; values[3*i+2] = diff;
			}

			rowPtr[_n-1] = count-3;
			colInd[count-3] = 0;    values[count-3] = diff;
			colInd[count-2] = _n-2; values[count-2] = diff;
			colInd[count-1] = _n-1; values[count-1] = -2*diff;

			jac = CSRMat<FP>(values, colInd, rowPtr, _n, _n, count);

			delete [] values;
			delete [] colInd;
			delete [] rowPtr;
		} else if( split == 2 ) {
			long count = 2*_n;
			FP* values = new FP[count];
			long* colInd = new long[count];
			long* rowPtr = new long[_n];

			rowPtr[0] = 0;
			colInd[0] = 1;    values[0] = adv;
			colInd[1] = _n-1; values[1] = -adv;

			for( long i = 1; i < _n-1; i++ ) {
				rowPtr[i] = 2*i;
				colInd[2*i]   = i-1; values[2*i]   = -adv;
				colInd[2*i+1] = i+1; values[2*i+1] = adv;
			}

			rowPtr[_n-1] = count-2;
			colInd[count-2] = 0;    values[count-2] = adv;
			colInd[count-1] = _n-2; values[count-1] = -adv;

			jac = CSRMat<FP>(values, colInd, rowPtr, _n, _n, count);

			delete [] values;
			delete [] colInd;
			delete [] rowPtr;
		} else {
			throw Exception() << "Analytic split " << split << " is not defined.";
		}
	}

	void SplitMat1Sparse(const FP t, const Vec<FP>& y, CSRMat<FP>& mat) {
		JacAnalyticSparse(1, t, y, mat);
	}

	void SplitMat2Sparse(const FP t, const Vec<FP>& y, CSRMat<FP>& mat) {
		JacAnalyticSparse(2, t, y, mat);
	}

	// These templates are necessary so the RHS can be called by both
	// standard calls and ADOL-C
	template<class T>
	inline void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _n; i++ )
			yp[i] = _diffusivity*y.PeriodicCentralDifference(i, 2, _fdorder, 1./_n);
	}

	template<class T>	
	inline void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _n; i++ )
			yp[i] = _advection*y.PeriodicCentralDifference(i, 1, _fdorder, 1./_n);
	}

public:
	AdvectionDiffusion1D(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		params["tf"].SetFP(0.1);

		_advection = GetDefaultFP(params,"adv",1./10);
		_diffusivity = GetDefaultFP(params,"diff",1.);
		_n = GetDefaultLong(params,"N",64);	
		
		_initialCondition.Resize(_n);
		for( long i = 0; i < _n; i++ )
			_initialCondition[i] = sin(2*M_PI*i/_n);
	}

	LINK_TWOSPLIT
	IVP_NAME("1D Advection-Diffusion")
};

#endif

