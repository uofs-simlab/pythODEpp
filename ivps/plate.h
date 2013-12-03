#ifndef PLATE_H
#define PLATE_H

class PLATE : public BaseIVP {
	FP _omega;
	FP _sigma;

	FP GetU(const Vec<FP>& y, long i, long j) {
		if( i < 1 ) return 0;
		if( i > 8 ) return 0;
		if( j < 1 ) return 0;
		if( j > 5 ) return 0;
		return y[8*(j-1)+i-1];
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		for( long i = 0; i < 40; i++ )
			yp[40+i] = y[40];

		for( int i = 1; i <= 8; i++ ) {
			for( int j = 1; j <= 5; j++ ) {
				FP plate = 20*GetU(y,i,j) - 8*(GetU(y,i+1,j)+GetU(y,i-1,j)+GetU(y,i,j+1)+GetU(y,i,j-1))
						   + 2*(GetU(y,i+1,j+1) + GetU(y,i-1,j-1) + GetU(y,i+1,j-1) + GetU(y,i-1,j+1))
						   + GetU(y,i+2,j) + GetU(y,i-2,j) + GetU(y,i,j+2) + GetU(y,i,j-2);
				plate /= 16./6561.;
				FP f;
				FP x = 2*i/9.;
				switch( j ) {
				case 2:
				case 4:
					f = 200*(exp(-5*sqr(t-x-2)) + exp(0-5*sqr(t-x-5)));
					break;
				default:
					f = 0;
				}
				yp[8*(j-1)+i-1] = -_omega*GetU(y, i, j) - _sigma*plate + f;
			}
		}
	
	}

public:
	PLATE(Hash<ParamValue>& params) : BaseIVP(params) {
		_omega = 1e3;
		_sigma = 1e2;

		params["tf"].SetFP(7.);

		// Initial conditions of the experiment
		_initialCondition.Resize(8*5*2);
		_initialCondition.Zero();
	}

	virtual const char* GetName() {
		return "PLATE";
	}
};

#endif

