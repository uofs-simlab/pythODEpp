#ifndef BRUSSELATOR_2D_H
#define BRUSSELATOR_2D_H

class Brusselator2D : public TwoSplittingIVP {
	long _n;
	FP _alpha;
	FP _B;
	FP _coeff;

	inline void JacValue(FP* values, long* colInd, long& pos, FP value, long ind) {
		values[pos] = value;
		colInd[pos] = ind;
		pos++;
	}

	inline void InnerOrder(long i, long* order) {
		if( i == 0 ) {
			order[0] = 2;
			order[1] = 3;
			order[2] = 4;
			order[3] = 1;
		} else if( i == _n-1 ) {
			order[0] = 4;
			order[1] = 1;
			order[2] = 2;
			order[3] = 3;
		} else {
			order[0] = 1;
			order[1] = 2;
			order[2] = 3;
			order[3] = 4;
		}
	}

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		if( split == 0 ) {
			long count = 12*_n*_n;
			FP* values = new FP[count];
			long* colInd = new long[count];
			long* rowPtr = new long[2*_n*_n];

			long pos = 0;
			for( long j = 0; j < _n; j++ ) {
				for( long i = 0; i < _n; i++ ) {
					long im1 = (i-1+_n)%_n;
					long ip1 = (i+1)%_n;
					long jm1 = (j-1+_n)%_n;
					long jp1 = (j+1)%_n;

					long order[6];
					if( j == 0 ) {
						InnerOrder(i, order);
						order[4] = 5;
						order[5] = 0;
					} else if( j == _n-1 ) {
						order[0] = 5;
						order[1] = 0;
						InnerOrder(i, order+2);
					} else {
						order[0] = 0;
						InnerOrder(i, order+1);
						order[5] = 5;
					}

					long posArray[] = {
						2*(_n*jm1+i), 2*(_n*j+im1), 2*(_n*j+i),
						2*(_n*j+i)+1, 2*(_n*j+ip1), 2*(_n*jp1+i)
					};

					FP valueArray[] = {
						_coeff, _coeff,
						-4*_coeff + 2*y[2*(j*_n+i)]*y[2*(j*_n+i)+1] - (_B+1),
						sqr(y[2*(j*_n+i)]),
						_coeff, _coeff };

					rowPtr[2*(j*_n+i)] = pos;
					for( long p = 0; p < 6; p++ )
						JacValue(values, colInd, pos, valueArray[order[p]], posArray[order[p]]);

					long posArray2[] = {
						2*(_n*jm1+i)+1, 2*(_n*j+im1)+1, 2*(_n*j+i),
						2*(_n*j+i)+1,   2*(_n*j+ip1)+1, 2*(_n*jp1+i)+1
					};
					valueArray[2] = _B - 2*y[2*(j*_n+i)]*y[2*(j*_n+i)+1];
					valueArray[3] = -4*_coeff - sqr(y[2*(j*_n+i)]);

					rowPtr[2*(j*_n+i)+1] = pos;
					for( long p = 0; p < 6; p++ )
						JacValue(values, colInd, pos, valueArray[order[p]], posArray2[order[p]]);
				}
			}

			jac = CSRMat<FP>(values, colInd, rowPtr, 2*_n*_n, 2*_n*_n, count);

			delete [] values;
			delete [] colInd;
			delete [] rowPtr;
		} else if( split == 1 ) {
			long count = 12*_n*_n;
			FP* values = new FP[count];
			long* colInd = new long[count];
			long* rowPtr = new long[2*_n*_n];

			long pos = 0;
			for( long j = 0; j < _n; j++ ) {
				for( long i = 0; i < _n; i++ ) {
					long im1 = (i-1+_n)%_n;
					long ip1 = (i+1)%_n;
					long jm1 = (j-1+_n)%_n;
					long jp1 = (j+1)%_n;

					long order[6];
					if( j == 0 ) {
						InnerOrder(i, order);
						order[4] = 5;
						order[5] = 0;
					} else if( j == _n-1 ) {
						order[0] = 5;
						order[1] = 0;
						InnerOrder(i, order+2);
					} else {
						order[0] = 0;
						InnerOrder(i, order+1);
						order[5] = 5;
					}

					long posArray[] = {
						2*(_n*jm1+i), 2*(_n*j+im1), 2*(_n*j+i),
						2*(_n*j+i)+1, 2*(_n*j+ip1), 2*(_n*jp1+i)
					};

					FP valueArray[] = {
						_coeff, _coeff,
						-4*_coeff, 0,
						_coeff, _coeff };

					rowPtr[2*(j*_n+i)] = pos;
					for( long p = 0; p < 6; p++ )
						JacValue(values, colInd, pos, valueArray[order[p]], posArray[order[p]]);

					long posArray2[] = {
						2*(_n*jm1+i)+1, 2*(_n*j+im1)+1, 2*(_n*j+i),
						2*(_n*j+i)+1,   2*(_n*j+ip1)+1, 2*(_n*jp1+i)+1
					};
					valueArray[2] = 0;
					valueArray[3] = -4*_coeff;

					rowPtr[2*(j*_n+i)+1] = pos;
					for( long p = 0; p < 6; p++ )
						JacValue(values, colInd, pos, valueArray[order[p]], posArray2[order[p]]);
				}
			}

			jac = CSRMat<FP>(values, colInd, rowPtr, 2*_n*_n, 2*_n*_n, count);

			delete [] values;
			delete [] colInd;
			delete [] rowPtr;
		} else {
			throw Exception() << "Analytic split " << split << " is not defined.";
		}
	}

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				long im1 = (i-1+_n)%_n;
				long ip1 = (i+1)%_n;
				long jm1 = (j-1+_n)%_n;
				long jp1 = (j+1)%_n;

				T uDxx = y[2*(j*_n+im1)] - 2*y[2*(j*_n+i)] + y[2*(j*_n+ip1)];
				T uDyy = y[2*(jm1*_n+i)] - 2*y[2*(j*_n+i)] + y[2*(jp1*_n+i)];
				T vDxx = y[2*(j*_n+im1)+1] - 2*y[2*(j*_n+i)+1] + y[2*(j*_n+ip1)+1];
				T vDyy = y[2*(jm1*_n+i)+1] - 2*y[2*(j*_n+i)+1] + y[2*(jp1*_n+i)+1];

				yp[2*(j*_n+i)]   = _coeff*(uDxx + uDyy);
				yp[2*(j*_n+i)+1] = _coeff*(vDxx + vDyy);
			}
		}
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				yp[2*(j*_n+i)]   = 1. + sqr(y[2*(j*_n+i)])*y[2*(j*_n+i)+1] - (_B+1)*y[2*(j*_n+i)];
				yp[2*(j*_n+i)+1] = _B*y[2*(j*_n+i)] - sqr(y[2*(j*_n+i)])*y[2*(j*_n+i)+1];

				FP x = FP(i)/_n;
				FP y = FP(j)/_n;
				if( sqr(x-0.3) + sqr(y-0.6) <= sqr(0.1) )
					yp[2*(j*_n+i)] += ProcessConditional(1.1-t, 0., 5.);
			}
		}
	}

public:
	Brusselator2D(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		SetDefaultFP(params, "tf", 11.5);

		_alpha = GetDefaultFP(params, "alpha", 1e-1);
		_n = GetDefaultLong(params, "N", 15);
		_coeff = _alpha*sqr(_n);
		_B = 3.4;
		
		_initialCondition.Resize(2*_n*_n);
		for( long i = 0; i < _n; i++ ) {
			for( long j = 0; j < _n; j++ ) {
				FP x = FP(i)/_n;
				FP y = FP(j)/_n;

				_initialCondition[2*(_n*j+i)]   = 22.*y*pow(1-y,3./2);
				_initialCondition[2*(_n*j+i)+1] = 27.*x*pow(1-x,3./2);
			}
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("2D Brusselator")
};

#endif

