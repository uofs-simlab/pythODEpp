#ifndef HEATTRANSFER_H
#define HEATTRANSFER_H

class HeatTransfer : public TwoSplittingIVP {
	long _nx;
	long _ny;
	FP _d;
	FP _uAvg;
	FP _length;
	FP _height;
	FP _tempIn;
	FP _wall;
	FP _dx;
	FP _dy;

	inline void JacValue(FP* values, long* colInd, long& pos, FP value, long ind) {
		values[pos] = value;
		colInd[pos] = ind;
		pos++;
	}

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
		if( split > 2 )
			throw Exception() << "Analytic split " << split << " is not defined.";

		long count = 5*_nx*_ny - 2*_nx - _ny;

		FP* values = new FP[count];
		long* colInd = new long[count];
		long* rowPtr = new long[_nx*_ny];

		long pos = 0;
		for( long j = 0; j < _ny; j++ ) {
			for( long i = 0; i < _nx; i++ ) {
				FP y = (j+0.5)*_dy - _height/2;
				FP v = -(3./2) * _uAvg*(1-sqr(y/(_height/2)))/(2*_dx);
				FP diff = _d;
				FP boundary = (j == 0 || j == _ny-1) ? 1 : 0;

				if( split == 1 ) v = 0;
				if( split == 2 ) diff = 0;

				rowPtr[j*_nx+i] = pos;
				if( j > 0 ) JacValue(values, colInd, pos, diff/sqr(_dy), (j-1)*_nx+i);
				if( i < _nx-1 ) {
					if( i > 0 ) JacValue(values, colInd, pos, diff/sqr(_dx) - v, j*_nx+i-1);
					JacValue(values, colInd, pos, -diff*(2/sqr(_dx)+(2-boundary)/sqr(_dy)) + v, j*_nx+i);
					JacValue(values, colInd, pos, diff/sqr(_dx), j*_nx+i+1);
				} else {
					JacValue(values, colInd, pos, -diff/sqr(_dx), j*_nx+i-2);
					JacValue(values, colInd, pos, 2*diff/sqr(_dx) - v, j*_nx+i-1);
					JacValue(values, colInd, pos, -diff*(1/sqr(_dx)+(2-boundary)/sqr(_dy)) + v, j*_nx+i);
				}
				if( j < _ny-1 ) JacValue(values, colInd, pos, diff/sqr(_dy), (j+1)*_nx+i);
			}
		}

		jac = CSRMat<FP>(values, colInd, rowPtr, _nx*_ny, _nx*_ny, count);

		delete [] values;
		delete [] colInd;
		delete [] rowPtr;
	}

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _nx; i++ ) {
			for( long j = 0; j < _ny; j++ ) {
				T ij   = y[_nx*j+i];
				T im1j = (i == 0)     ? _tempIn : y[_nx*j+i-1];
				T ijm1 = (j == 0)     ? (_wall*_dx+ij) : y[_nx*(j-1)+i];
				T ijp1 = (j == _ny-1) ? (_wall*_dx+ij) : y[_nx*(j+1)+i];

				if( i == _nx-1 )
					yp[_nx*j+i] = _d*((-y[_nx*j+i-2]+2*im1j-ij)/sqr(_dx) + (ijm1 + ijp1 - 2*ij)/sqr(_dy));
				else
					yp[_nx*j+i] = _d*((im1j+y[_nx*j+i+1]-2*ij)/sqr(_dx) + (ijm1 + ijp1 - 2*ij)/sqr(_dy));
			}
		}
	}
	
	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _nx; i++ ) {
			for( long j = 0; j < _ny; j++ ) {
				T ij   = y[_nx*j+i];
				T im1j = (i == 0) ? _tempIn : y[_nx*j+i-1];

				FP y = (j+0.5)*_dy - _height/2;
				yp[_nx*j+i] = -(3./2) * _uAvg*(1-sqr(y/(_height/2))) * (ij-im1j)/(2*_dx);
			}
		}
	}

public:
	HeatTransfer(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		SetDefaultFP(params, "tf", 50.);
	
		_nx     = GetDefaultLong(params, "NX", 200);
		_ny     = GetDefaultLong(params, "NY", 40);
		_d      = GetDefaultFP(params, "d", 1.38e-7);
		_length = GetDefaultFP(params, "height", 0.01);
		_height = GetDefaultFP(params, "length", 0.7);
		_tempIn = GetDefaultFP(params, "tempin", 30+273.15);
		_wall   = GetDefaultFP(params, "wall", 10000./0.58);
		_uAvg   = GetDefaultFP(params, "uavg", 1e-6*4.806*sqr(_height)/(12*8.9e-4));

		_dx = _length/_nx;
		_dy = _height/_ny;
		
		_initialCondition.Resize(_nx*_ny);
		for( long i = 0; i < _nx*_ny; i++ )
			_initialCondition[i] = 25+273.15;
	}

	LINK_TWOSPLIT
	IVP_NAME("Heat Transfer Equation");
};

#endif

