#ifndef RKP_H
#define RKP_H

class RKP : public TwoSplittingIVP {
	long _nx, _ny;
	FP _dx, _dy;
	FP _l1, _l2;
	FP _a, _b;
	FP _mu;
	
	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _nx; i++ ) {
			for( long j = 0; j < _ny; j++ ) {
				if( i == 0 )
					yp[j*_nx+i] = _l1*(-y[j*_nx+i] + y[j*_nx+i+1]) / (_dx*_dx);
				else if( i == _nx-1 )
					yp[j*_nx+i] = _l1*(y[j*_nx+i-1] - y[j*_nx+i]) / (_dx*_dx);
				else
					yp[j*_nx+i] = _l1*(y[j*_nx+i-1] - 2*y[j*_nx+i] + y[j*_nx+i+1]) / (_dx*_dx);

				if( j == 0 )
					yp[j*_nx+i] += _l2*(-y[j*_nx+i] + y[(j+1)*_nx+i]) / (_dy*_dy);
				else if( j == _ny-1 )
					yp[j*_nx+i] += _l2*(y[(j-1)*_nx+i] - y[j*_nx+i]) / (_dy*_dy);
				else
					yp[j*_nx+i] += _l2*(y[(j-1)*_nx+i] - 2*y[j*_nx+i] + y[(j+1)*_nx+i]) / (_dy*_dy);
			}
		}
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp) {
		for( long i = 0; i < _nx; i++ ) {
			for( long j = 0; j < _ny; j++ ) {
				yp[j*_nx + i] = _mu*y[j*_nx+i]*(1-y[j*_nx+i]);
			}
		}
	}

public:
	RKP(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		params["tf"].SetFP(5.);
		
		_a = 2;
		_b = 1;
		_nx = 20;
		_ny = 20;
		_mu = 0.1;
		_l1 = 0.1;
		_l2 = 0.01;

		ParamValue* pv;
		if( (pv = params.Get("NX")) )
			_nx = pv->GetLong();
		if( (pv = params.Get("NY")) )
			_ny = pv->GetLong();
		if( (pv = params.Get("a")) )
			_a = pv->GetFP();
		if( (pv = params.Get("b")) )
			_b = pv->GetFP();
		if( (pv = params.Get("l1")) )
			_l1 = pv->GetFP();
		if( (pv = params.Get("l2")) )
			_l2 = pv->GetFP();
		if( (pv = params.Get("mu")) )
			_mu = pv->GetFP();

		_dx = 2*_a/(_nx-1);
		_dy = 2*_b/(_ny-1);

		_initialCondition.Resize(_nx*_ny);
		for( long i = 0; i < _nx; i++ ) {
			for( long j = 0; j < _ny; j++ ) {
				FP x = -_a + (i-1)*_dx;
				FP y = -_b + (i-1)*_dy;
				if( x*x + 4*y*y <= 0.25 )
					_initialCondition[j*_nx + i] = 1.;
				else
					_initialCondition[j*_nx + i] = exp(-10*(x*x + 4*y*y - 0.25));
			}
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("RKP Equation")
};

#endif

