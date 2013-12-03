#ifndef SCOTT_WANG_SHOWALTER_H
#define SCOTT_WANG_SHOWALTER_H

class ScottWangShowalter : public TwoSplittingIVP {
protected:
	long _n;
	FP _size;
	FP _h;
	FP _mu;
	FP _Le;
	FP _eps;
	FP _kappa;

	enum Problem {
		Problem1,
		Problem2,
		Problem3
	};
	Problem _problem;
	
	inline void JacValue(FP* values, long* colInd, long& pos, FP value, long ind) {
		values[pos] = value;
		colInd[pos] = ind;
		pos++;
	}

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		if( split > 1 ) 
			throw Exception() << "Analytic split " << split << " is not defined.";

		long count = 12*_n*_n - 8*_n;
		FP* values = new FP[count];
		long* colInd = new long[count];
		long* rowPtr = new long[2*_n*_n];

		long pos = 0;
		for( long j = 0; j < _n; j++ ) {
			for( long i = 0; i < _n; i++ ) {
				FP Aij = y[2*(j*_n+i)];
				FP Tij = y[2*(j*_n+i)+1];
				FP f = exp(Tij/(1+_eps*Tij));

				FP advij   = 0;
				FP advijm1 = 0, advijp1 = 0;
				FP advim1j = 0, advip1j = 0;

				if( split == 0 && _problem != Problem1 ) {
					FP v1, v2;
					if( _problem == Problem2 ) {
						v1 = 80;
						v2 = 80;
					} else {
						FP x = (i+1)*_h;
						FP y = (j+1)*_h;
						v1 = -50*(y-5);
						v2 =  50*(x-5);
					}

					v1 /= 2*_h;
					v2 /= 2*_h;

					if( v1 > 0 ) {
						advij   = v1;
						advim1j = i == 0 ? 0 : -v1;
					} else {
						advip1j = i == _n-1 ? 0 : v1;
						advij   = -v1;
					}

					if( v2 > 0 ) {
						advij  += v2;
						advijm1 = j == 0 ? 0 : -v2;
					} else {
						advijp1 = j == _n-1 ? 0 : v2;
						advij  += -v2;
					}
				}

				// First row
				rowPtr[2*(j*_n+i)] = pos;
				if( j > 0 )
					JacValue(values, colInd, pos, 1/sqr(_h) - advijm1, 2*((j-1)*_n+i));
				if( i > 0 )
					JacValue(values, colInd, pos, 1/sqr(_h) - advim1j, 2*(j*_n+i-1));
				if( split == 0 ) {
					JacValue(values, colInd, pos, -4/sqr(_h) - f - advij, 2*(j*_n+i));
					JacValue(values, colInd, pos, -Aij*f/sqr(_eps*Tij+1), 2*(j*_n+i)+1);
				} else {
					JacValue(values, colInd, pos, -4/sqr(_h), 2*(j*_n+i));
					JacValue(values, colInd, pos, 0, 2*(j*_n+i)+1);
				}
				if( i < _n-1 )
					JacValue(values, colInd, pos, 1/sqr(_h) - advip1j, 2*(j*_n+i+1));
				if( j < _n-1 )
					JacValue(values, colInd, pos, 1/sqr(_h) - advijp1, 2*((j+1)*_n+i));

				// Second row
				rowPtr[2*(j*_n+i)+1] = pos;
				if( j > 0 )
					JacValue(values, colInd, pos, _Le/sqr(_h) - advijm1, 2*((j-1)*_n+i)+1);
				if( i > 0 )
					JacValue(values, colInd, pos, _Le/sqr(_h) - advim1j, 2*(j*_n+i-1)+1);
				if( split == 0 ) {
					JacValue(values, colInd, pos, (1/_kappa)*f, 2*(j*_n+i));
					JacValue(values, colInd, pos, -4*_Le/sqr(_h) + (1/_kappa)*(Aij*f/sqr(_eps*Tij+1) - 1) - advij, 2*(j*_n+i)+1);
				} else {
					JacValue(values, colInd, pos, 0, 2*(j*_n+i));
					JacValue(values, colInd, pos, -4*_Le/sqr(_h), 2*(j*_n+i)+1);
				}
				if( i < _n-1 )
					JacValue(values, colInd, pos, _Le/sqr(_h) - advip1j, 2*(j*_n+i+1)+1);
				if( j < _n-1 )
					JacValue(values, colInd, pos, _Le/sqr(_h) - advijp1, 2*((j+1)*_n+i)+1);
			}
		}

		jac = CSRMat<FP>(values, colInd, rowPtr, 2*_n*_n, 2*_n*_n, count);

		delete [] values;
		delete [] colInd;
		delete [] rowPtr;
	}

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		FP Ab = _mu / exp(_mu/(1+_eps*_mu));
		FP Tb = _mu;

		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				T Aij   = y[2*(j*_n+i)];
				T Aim1j = i == 0    ? Ab : y[2*(j*_n+i-1)];
				T Aip1j = i == _n-1 ? Ab : y[2*(j*_n+i+1)];
				T Aijm1 = j == 0    ? Ab : y[2*((j-1)*_n+i)];
				T Aijp1 = j == _n-1 ? Ab : y[2*((j+1)*_n+i)];
				
				T Tij   = y[2*(j*_n+i)+1];
				T Tim1j = i == 0    ? Tb : y[2*(j*_n+i-1)+1];
				T Tip1j = i == _n-1 ? Tb : y[2*(j*_n+i+1)+1];
				T Tijm1 = j == 0    ? Tb : y[2*((j-1)*_n+i)+1];
				T Tijp1 = j == _n-1 ? Tb : y[2*((j+1)*_n+i)+1];

				T lapA = (Aim1j + Aip1j + Aijm1 + Aijp1 - 4*Aij) / sqr(_h);
				T lapT = (Tim1j + Tip1j + Tijm1 + Tijp1 - 4*Tij) / sqr(_h);

				yp[2*(j*_n+i)]   = lapA;
				yp[2*(j*_n+i)+1] = _Le*lapT;
			}
		}
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		FP Ab = _mu / exp(_mu/(1+_eps*_mu));
		FP Tb = _mu;

		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				T Aij   = y[2*(j*_n+i)];
				T Tij   = y[2*(j*_n+i)+1];
				T f = exp(Tij/(1+_eps*Tij));

				yp[2*(j*_n+i)]   = _mu - Aij*f;
				yp[2*(j*_n+i)+1] = (1/_kappa)*(Aij*f-Tij);

				if( _problem != Problem1 ) {
					FP v1, v2;
					if( _problem == Problem2 ) {
						v1 = 80;
						v2 = 80;
					} else {
						FP x = (i+1)*_h;
						FP y = (j+1)*_h;
						v1 = -50*(y-5);
						v2 =  50*(x-5);
					}

					T Adx, Ady;
					T Tdx, Tdy;

					if( v1 > 0 ) {
						T Aim1j = i == 0 ? Ab : y[2*(j*_n+i-1)];
						T Tim1j = i == 0 ? Tb : y[2*(j*_n+i-1)+1];
						Adx = (Aij-Aim1j)/(2*_h);
						Tdx = (Tij-Tim1j)/(2*_h);
					} else {
						T Aip1j = i == _n-1 ? Ab : y[2*(j*_n+i+1)];
						T Tip1j = i == _n-1 ? Tb : y[2*(j*_n+i+1)+1];
						Adx = (Aip1j-Aij)/(2*_h);
						Tdx = (Tip1j-Tij)/(2*_h);
					}
					
					if( v2 > 0 ) {
						T Aijm1 = j == 0 ? Ab : y[2*((j-1)*_n+i)];
						T Tijm1 = j == 0 ? Tb : y[2*((j-1)*_n+i)+1];
						Ady = (Aij-Aijm1)/(2*_h);
						Tdy = (Tij-Tijm1)/(2*_h);
					} else {
						T Aijp1 = j == _n-1 ? Ab : y[2*((j+1)*_n+i)];
						T Tijp1 = j == _n-1 ? Tb : y[2*((j+1)*_n+i)+1];
						Ady = (Aijp1-Aij)/(2*_h);
						Tdy = (Tijp1-Tij)/(2*_h);
					}

					yp[2*(j*_n+i)]   -= v1*Adx + v2*Ady;
					yp[2*(j*_n+i)+1] -= v1*Tdx + v2*Tdy;
				}
			}
		}
	}
 
public:
	ScottWangShowalter(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		SetDefaultFP(params, "tf", 0.06);
		
		_n = GetDefaultLong(params, "N", 20);

		_problem = Problem1;
		ParamValue* pv;
		if( (pv = params.Get("problem")) ) {
			if( std::string(pv->GetString()) == "Problem1" )
				_problem = Problem1;
			else if( std::string(pv->GetString()) == "Problem2" )
				_problem = Problem2;
			else if( std::string(pv->GetString()) == "Problem3" )
				_problem = Problem3;
		}

		_size = 10.;
		_h = _size/(_n+1);
		_mu = 1.8;
		_Le = 1.;
		_eps = 0.18;
		_kappa = 0.0005;

		_initialCondition.Resize(2*_n*_n);
		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				FP x = (i+1)*_h;
				FP y = (j+1)*_h;

				_initialCondition[2*(j*_n+i)]   = _mu / exp(_mu/(1+_eps*_mu));

				switch( _problem ) {
				case Problem1:
					_initialCondition[2*(j*_n+i)+1] = _mu + exp(-50*(sqr(x-5) + sqr(y-(y<x?3.5:6.5))));
					break;
				case Problem2:
					_initialCondition[2*(j*_n+i)+1] = _mu + exp(-50*(sqr(x-(y<x?4:3.5)) + sqr(y-(y<x?3.5:4))));
				case Problem3:
					_initialCondition[2*(j*_n+i)+1] = _mu + exp(-50*(sqr(x-5) + sqr(y-(y<x?3:7))));
					break;
				}
			}
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("Scott-Wang-Showalter");
};

#endif

