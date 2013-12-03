#ifndef COMBUSTION_H
#define COMBUSTION_H

class CombustionARD : public TwoSplittingIVP {
	long _n, _m;
	FP _x0, _xmin, _xmax, _dx;
	FP _L, _U0;
	FP _gamma, _beta, _sigma;
	FP _alpha1, _alpha2, _alpha3;
	FP _cs;

	enum Reaction {
		FKPP,
		Ignition,
		Fisher
	};

	Reaction _reaction;

	inline void JacValue(FP* values, long* colInd, long& pos, FP value, long ind) {
		values[pos] = value;
		colInd[pos] = ind;
		pos++;
	}

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		if( split > 2 ) 
			throw Exception() << "Analytic split " << split << " is not defined.";

		long count = 3*_n;
		FP* values = new FP[count];
		long* colInd = new long[count];
		long* rowPtr = new long[_n];

		long pos = 0;
		for( long i = 0; i < _n; i++ ) {
			FP x = _xmin + i*_dx;

			// Diffusion
			FP d = split != 2 ? (1+_U0*sin(M_PI*x/_L))*(_gamma/_beta)/sqr(_dx) : 0;

			// Advection
			FP a = split != 1 ? -(1+_U0*sin(M_PI*x/_L))/(2*_dx) : 0;

			// Reaction
			FP r = 0;
			if( split != 1 ) {
				switch( _reaction ) {
				case FKPP:
					r = _alpha1*(1-2*y[i]);
					break;
				case Ignition:
					r = _cs > y[i] ? 0 : _alpha2*(_cs-2*y[i]+1);
					break;
				case Fisher:
					r = -_alpha3*pow(y[i], _m-1.)*(_m*(y[i]-1)+y[i]);
					break;
				}
			}

			rowPtr[i] = pos;
			if( i == 0 ) {
				JacValue(values, colInd, pos, -2*d+r, i);
				JacValue(values, colInd, pos, d+a,    i+1);
				JacValue(values, colInd, pos, d-a,    _n-1);
			} else if( i == _n-1 ) {
				JacValue(values, colInd, pos, d+a,    0);
				JacValue(values, colInd, pos, d-a,    i-1);
				JacValue(values, colInd, pos, -2*d+r, i);
			} else {
				JacValue(values, colInd, pos, d-a,    i-1);
				JacValue(values, colInd, pos, -2*d+r, i);
				JacValue(values, colInd, pos, d+a,    i+1);
			}
		}

		jac = CSRMat<FP>(values, colInd, rowPtr, _n, _n, count);

		delete [] values;
		delete [] colInd;
		delete [] rowPtr;
	}

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _n; i++ ) {
			FP x = _xmin + i*_dx;
			yp[i] = (1+_U0*sin(M_PI*x/_L))*(_gamma/_beta)*y.PeriodicCentralDifference(i, 2, 2, _dx);
		}
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _n; i++ ) {
			FP x = _xmin + i*_dx;

			switch( _reaction ) {
			case FKPP:
				yp[i] = _alpha1*y[i]*(1-y[i]);
				break;
			case Ignition:
				yp[i] = ProcessConditional(_cs - y[i], 0, _alpha2*(y[i]-_cs)*(1-y[i]));
				break;
			case Fisher:
				yp[i] = _alpha3*pow(y[i], (FP)_m)*(1-y[i]);
				break;
			}

			yp[i] -= (1+_U0*sin(M_PI*x/_L))*y.PeriodicCentralDifference(i, 1, 2, _dx);
		}
	}

public:
	CombustionARD(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		SetDefaultFP(params, "tf", 30.);
				
		ParamValue* pv;
		_reaction = FKPP;
		if( (pv = params.Get("reaction")) ) {
			if( std::string("FKPP") == pv->GetString() )
				_reaction = FKPP;
			else if( std::string("Ignition") == pv->GetString() )
				_reaction = Ignition;
			else if( std::string("Fisher") == pv->GetString() )
				_reaction = Fisher;
		}

		_n = 40;
		if( (pv = params.Get("N")) )
			_n = pv->GetLong();
		_U0 = 0.75;
		if( (pv = params.Get("U0")) )
			_U0 = pv->GetFP();	

		_m = 10;
		_x0 = 20.5;
		_xmin = 10.;
		_xmax = 50.;
		_dx = (_xmax-_xmin)/_n;
		_beta = 1;
		_cs = 0.6;
		_gamma = 0.1;
		_alpha1 = 0.1;
		_alpha2 = _alpha1/sqr(1-_cs);
		_alpha3 = _alpha1/(4*(pow(_m/(_m+1.), (FP)_m) * (1 - _m/(_m+1.))));
		_sigma = 10.;
		_L = 1;

		_initialCondition.Resize(_n);
		for( long i = 0; i < _n; i++ ) {
			FP x = _xmin + i*_dx;
			_initialCondition[i] = exp(-sqr(x-_x0)/_sigma);
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("Combustion ARD")
};

#endif

