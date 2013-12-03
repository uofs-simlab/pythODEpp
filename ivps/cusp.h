#ifndef CUSP_H
#define CUSP_H

class CUSP : public TwoSplittingIVP {
	long _N;
	FP _sigma;
	FP _eps;

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _N; i++ ) {
			long ip1 = (i+1)%_N;
			long im1 = (_N+i-1)%_N;
			T D = _sigma*_N*_N;

			yp[3*i]   = D*(y[3*im1]   - 2*y[3*i]   + y[3*ip1]);
			yp[3*i+1] = D*(y[3*im1+1] - 2*y[3*i+1] + y[3*ip1+1]);
			yp[3*i+2] = D*(y[3*im1+2] - 2*y[3*i+2] + y[3*ip1+2]);
		}
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _N; i++ ) {
			T yi = y[3*i];
			T ai = y[3*i+1];
			T bi = y[3*i+2];
			T u  = (yi-0.7)*(yi-1.3);
			T v  = u/(u+0.1);

			yp[3*i]   = -(1./_eps)*(yi*yi*yi + ai*yi + bi);
			yp[3*i+1] = bi + 0.07*v;
			yp[3*i+2] = (1-sqr(ai))*bi - ai - 0.4*yi + 0.035*v;
		}
	//	yp.Zero();
	}

public:
	CUSP(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		_N = GetDefaultLong(params,"N",32);
		_sigma = GetDefaultFP(params, "sigma", 1./144);
		_eps = GetDefaultFP(params,"eps",1e-4);

		SetDefaultFP(params,"tf",1.1);

		// Initial conditions of the experiment
		_initialCondition.Resize(3*_N);
		for( int i = 0; i < _N; i++ ) {
			_initialCondition[3*i]   = 0.;
			_initialCondition[3*i+1] = -2*cos(2.*i*M_PI/_N);
			_initialCondition[3*i+2] = 2*sin(2.*i*M_PI/_N);
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("CUSP")
};

#endif

