#ifndef CHEMOTAXIS_H
#define CHEMOTAXIS_H

class Angiogenesis1D : public TwoSplittingIVP {
	long _n;
	FP _eps;
	FP _delta;
	FP _alpha;
	FP _beta;
	FP _gamma;
	FP _kappa;
	FP _lambda;
	FP _mu;
	FP _cstar;
	FP _dx;

	inline void JacValue(FP* values, long* colInd, long& pos, FP value, long ind) {
		values[pos] = value;
		colInd[pos] = ind;
		pos++;
	}

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		if( split > 2 ) 
			throw Exception() << "Analytic split " << split << " is not defined.";

		long count = 10*_n - 6;
		FP* values = new FP[count];
		long* colInd = new long[count];
		long* rowPtr = new long[2*_n];

		FP e = split == 2 ? 0 : _eps;
		FP d = split == 2 ? 0 : _delta;
		FP k = split == 1 ? 0 : _kappa;
		FP r = split == 1 ? 0 : 1;

		long pos = 0;
		for( long i = 0; i < _n; i++ ) {
			FP p   = y[2*i];
			FP c   = y[2*i+1];
			FP pm1 = i > 0 ? y[2*(i-1)] : 0;
			FP pp1 = i < _n-1 ? y[2*(i+1)] : 1;
			FP cm1 = i > 0 ? y[2*(i-1)+1] : 1;
			FP cp1 = i < _n-1 ? y[2*(i+1)+1] : 0;

			// First row
			rowPtr[2*i] = pos;
			if( i > 0 ) {
				//FP dpidpim1 = k*(cp1-cm1)/sqr(2*_dx);
				//FP dpicpim1 = -k*(p/sqr(_dx) - (pp1-pm1)/sqr(2*_dx));
				FP dpidpim1 = (-k/(2*sqr(_dx)))*(cm1-c);
				FP dpicpim1 = (-k/(2*sqr(_dx)))*(pm1+p);
				JacValue(values, colInd, pos, e/sqr(_dx) + dpidpim1, 2*(i-1));
				JacValue(values, colInd, pos, dpicpim1, 2*(i-1)+1);
			}
			//FP dpidpi = -k*(cm1-2*c+cp1)/sqr(_dx);
			//FP dpicpi = 2*k*p/sqr(_dx);
			FP dpidpi = (-k/(2*sqr(_dx)))*(cm1-2*c+cp1);
			FP dpicpi = (-k/(2*sqr(_dx)))*(-1)*(pm1+2*p+pp1);
			JacValue(values, colInd, pos, -2*e/sqr(_dx) + dpidpi + r*(_mu*(1-2*p)*fmax(0.,c-_cstar) - _beta), 2*i);
			JacValue(values, colInd, pos, dpicpi + r*(_mu*p*(1-p)*(c-_cstar>0?1:0)), 2*i+1);
			if( i < _n-1 ) {
				//FP dpidpip1 = -k*(cp1-cm1)/sqr(2*_dx);
				//FP dpicpip1 = -k*(p/sqr(_dx) + (pp1-pm1)/sqr(2*_dx));
				FP dpidpip1 = (-k/(2*sqr(_dx)))*(cp1-c);
				FP dpicpip1 = (-k/(2*sqr(_dx)))*(p+pp1);
				JacValue(values, colInd, pos, e/sqr(_dx) + dpidpip1, 2*(i+1));
				JacValue(values, colInd, pos, dpicpip1, 2*(i+1)+1);
			}

			// Second row
			rowPtr[2*i+1] = pos;
			if( i > 0 ) JacValue(values, colInd, pos, d/sqr(_dx), 2*(i-1)+1);
			JacValue(values, colInd, pos, r*(-_alpha*c/(_gamma+c)), 2*i);
			JacValue(values, colInd, pos, -2*d/sqr(_dx) - r*(_lambda + _alpha*p*_gamma/sqr(_gamma+c)), 2*i+1);
			if( i < _n-1 ) JacValue(values, colInd, pos, d/sqr(_dx), 2*(i+1)+1);
		}

		jac = CSRMat<FP>(values, colInd, rowPtr, 2*_n, 2*_n, count);

		delete [] values;
		delete [] colInd;
		delete [] rowPtr;
	}

	template <class T>	
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		T pxx0 = (-2*y[0] + y[2])/sqr(_dx);
		T cxx0 = (1 - 2*y[1] + y[3])/sqr(_dx);
		yp[0] = _eps*pxx0;
		yp[1] = _delta*cxx0;

		for( int i = 1; i < _n-1; i++ ) {
			T pxx = y.PeriodicCentralDifference(2*i,   2, 2, _dx, 2);
			T cxx = y.PeriodicCentralDifference(2*i+1, 2, 2, _dx, 2);
			yp[2*i]   = _eps*pxx;
			yp[2*i+1] = _delta*cxx;
		}

		T pxxn = (y[2*(_n-2)]   - 2*y[2*(_n-1)] + 1)/sqr(_dx);
		T cxxn = (y[2*(_n-2)+1] - 2*y[2*(_n-1)+1])/sqr(_dx);
		yp[2*(_n-1)]   = _eps*pxxn;
		yp[2*(_n-1)+1] = _delta*cxxn;
	}

	template <class T>	
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		// Advection
		/*T cxx0 = (1 - 2*y[1] + y[3])/sqr(_dx);
		T px0  = y[2]/(2*_dx);
		T cx0  = (y[3]-1)/(2*_dx);
		yp[0] = -_kappa*(cxx0*y[0] + cx0*px0);
		yp[1] = 0;

		for( int i = 1; i < _n-1; i++ ) {
			T cxx = y.PeriodicCentralDifference(2*i+1, 2, 2, _dx, 2);
			T px  = y.PeriodicCentralDifference(2*i,   1, 2, _dx, 2);
			T cx  = y.PeriodicCentralDifference(2*i+1, 1, 2, _dx, 2);
			yp[2*i]   = -_kappa*(cxx*y[2*i] + cx*px);
			yp[2*i+1] = 0;
		}

		T cxxn = (y[2*(_n-2)+1] - 2*y[2*(_n-1)+1])/sqr(_dx);
		T pxn  = (1 - y[2*(_n-2)])/(2*_dx); 
		T cxn  = -y[2*(_n-2)+1]/(2*_dx);
		yp[2*(_n-1)]   = -_kappa*(cxxn*y[2*(_n-1)] + cxn*pxn);
		yp[2*(_n-1)+1] = 0;*/
		for( long i = 0; i < _n; i++ ) {
			T pi = y[2*i];
			T pim1 = (i==0)    ? 0 : y[2*(i-1)];
			T pip1 = (i==_n-1) ? 1 : y[2*(i+1)];
			T ci = y[2*i+1];
			T cim1 = (i==0)    ? 1 : y[2*(i-1)+1];
			T cip1 = (i==_n-1) ? 0 : y[2*(i+1)+1];

			yp[2*i]   = -_kappa/(2*sqr(_dx)) * (cip1*(pi+pip1) - ci*(pim1+2*pi+pip1) + cim1*(pim1+pi));
			yp[2*i+1] = 0;
		}

		// Reaction
		for( int i = 0; i < _n; i++ ) {
			T p = y[2*i];
			T c = y[2*i+1];
			T diff = c - _cstar;

			yp[2*i]   += _mu*p*(1-p)*CustomMax0(diff) - _beta*p;
			yp[2*i+1] += -_lambda*c - _alpha*p*c/(_gamma+c);
		}
	}

public:
	Angiogenesis1D(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		SetDefaultFP(params, "tf", 0.7);

		_n = 100;
		_eps = 1e-3;
		_delta = 1.;
		_alpha = 10.;
		_beta = 4.;
		_gamma = 1.;
		_kappa = 0.75;
		_lambda = 1.;
		_mu = 100.;
		_cstar = 0.2;

		ParamValue* pv;		
		if( (pv = params.Get("N")) ) _n = pv->GetLong();
		if( (pv = params.Get("d")) ) _delta = pv->GetFP();

		_dx = 1./(_n+1);

		_initialCondition.Resize(2*_n);
		for( long i = 0; i < _n; i ++ ) {
			_initialCondition[2*i]   = 0.;
			_initialCondition[2*i+1] = cos(0.5*M_PI*_dx*(i+1));
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("Tumour Angiogenesis 1D")
};

#endif

