#ifndef BURGERS_H
#define BURGERS_H

class Burgers : public TwoSplittingIVP {
	long _n;
	FP _h;
	FP _d;

	inline void JacValue(FP* values, long* colInd, long& pos, FP value, long ind) {
		values[pos] = value;
		colInd[pos] = ind;
		pos++;
	}

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		if( split > 2 ) 
			throw Exception() << "Analytic split " << split << " is not defined.";

		long count = 5*_n - 8;
		FP* values = new FP[count];
		long* colInd = new long[count];
		long* rowPtr = new long[_n];

		FP c0 = -5./2;
		FP c1 =  4./3;
		FP c2 = -1./12;
		FP d = split != 2 ? _d/(_h*_h) : 0;
		FP a = split != 1 ? -1./_h : 0;

		long pos = 0;
		for( long i = 0; i < _n; i++ ) {
			FP yim1 = i == 0 ? ExactSolution<FP>(0,t) : y[i-1];
			FP yip1 = i == _n-1 ? ExactSolution<FP>(1,t) : y[i+1];
			FP yi   = y[i];

			FP am1 = y[i] > 0 ?  -yi : 0;
			FP ap1 = y[i] > 0 ?   0    : yi;
			FP ai  = y[i] > 0 ? (2*yi-yim1) : (yip1-2*yi);

			rowPtr[i] = pos;
			if( i == 0 ) {
				JacValue(values, colInd, pos, -2*d + ai*a, i);
				JacValue(values, colInd, pos,  d + ap1*a,   i+1);
			} else if( i == _n-1 ) {
				JacValue(values, colInd, pos,  d + am1*a,   i-1);
				JacValue(values, colInd, pos, -2*d + ai*a, i);
			} else {
				if( i > 1 ) JacValue(values, colInd, pos, c2*d, i-2);
				JacValue(values, colInd, pos, c1*d + am1*a, i-1);
				JacValue(values, colInd, pos, c0*d + ai*a,  i);
				JacValue(values, colInd, pos, c1*d + ap1*a, i+1);
				if( i < _n-2 ) JacValue(values, colInd, pos, c2*d, i+2);
			}
		}

		jac = CSRMat<FP>(values, colInd, rowPtr, _n, _n, count);

		delete [] values;
		delete [] colInd;
		delete [] rowPtr;
	}

	template <class T>
	T ExactSolution(T x, T t) {
		T r1 = exp(-(20*(x-0.5)+99*t)/(400*_d));
		T r2 = exp(-(4*(x-0.5)+3*t)/(16*_d));
		T r3 = exp(-(x-3./8)/(2*_d));
		return 1. - (0.9*r1+0.5*r2) / (r1+r2+r3);
	}

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		FP coeff = _d/(_h*_h);
		FP c0 = -5./2;
		FP c1 =  4./3;
		FP c2 = -1./12;

		T bc0 = ExactSolution<T>(0,t);
		T bc1 = ExactSolution<T>(1,t);

		yp[0] = coeff*(bc0 - 2*y[0] + y[1]);
		yp[1] = coeff*(c2*bc0 + c1*y[0] + c0*y[1] + c1*y[2] + c2*y[3]);

		for( long i = 2; i < _n-2; i++ )
			yp[i] = coeff*(c2*y[i-2] + c1*y[i-1] + c0*y[i] + c1*y[i+1] + c2*y[i+2]);
		
		yp[_n-2] = coeff*(c2*y[_n-4] + c1*y[_n-3] + c0*y[_n-2] + c1*y[_n-1] + c2*bc1);
		yp[_n-1] = coeff*(y[_n-2] - 2*y[_n-1] + bc1);
	}

	inline FP UpwindConditional(FP yim1, FP yi, FP yip1) {
		return (yi > 0) ? (yi - yim1) : (yip1 - yi);
	}

#ifdef USE_ADOL_C
	inline adouble UpwindConditional(adouble yim1, adouble yi, adouble yip1) {
		adouble ret;
		condassign(ret, yi, yi-yim1, yip1-yi);
		return ret;
	}
#endif

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		FP coeff = 1/_h;

		T bc0 = ExactSolution<T>(0,t);
		T bc1 = ExactSolution<T>(1,t);

		yp[0] = -y[0]*coeff*UpwindConditional(bc0, y[0], y[1]);
		for( long i = 1; i < _n-1; i++ )
			yp[i] = -y[i]*coeff*UpwindConditional(y[i-1], y[i], y[i+1]);
		yp[_n-1] = -y[_n-1]*coeff*UpwindConditional(y[_n-2], y[_n-1], bc1);
	}

public:
	Burgers(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		SetDefaultFP(params, "tf", 3.);
		_n = GetDefaultLong(params, "N", 100);
		_h = 1./(_n+1);
		_d = GetDefaultFP(params, "d", 0.1);

		_initialCondition.Resize(_n);
		for( long i = 0; i < _n; i++ )
			_initialCondition[i] = ExactSolution<FP>((i+1)*_h, 0);
	}

	LINK_TWOSPLIT
	IVP_NAME("Burgers' Equation")
};

#endif

