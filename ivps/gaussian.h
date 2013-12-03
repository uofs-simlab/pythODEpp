#ifndef GAUSSIAN_H
#define GAUSSIAN_H

class Gaussian : public TwoSplittingIVP {
	long _n;
	FP _h;

public:
	Gaussian(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		params["tf"].SetFP(5.);
		_n = 50;

		ParamValue* pv;
		if( (pv = params.Get("N")) )
			_n = pv->GetLong();

		_h = 20./(_n+1);

		_initialCondition.Resize(_n);
		for( long i = 0; i < _n; i++ )
			_initialCondition[i] = exp(-sqr(-10 + (i+1)*_h));
	}

	void Split1(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		for( long i = 0; i < _n; i++ ) {
			//FP x = -10 + (i+1)*_h;
			FP yi   = y[i];
			FP yip1 = (i == 0)    ? exp(-10*10 - t) : y[i+1];
			FP yim1 = (i == _n-1) ? exp(-10*10 - t) : y[i-1];

			yp[i] = (1/(2*_h*_h)) * ((yip1+yi)*(yip1-yi) - (yim1+yi)*(yi-yim1));
		}
	}

	void Split2(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		for( long i = 0; i < _n; i++ ) {
			FP x = -10 + (i+1)*_h;
			yp[i] = -y[i] + (2-8*x*x)*exp(-2*(x*x+t));
		}
	}

	virtual const char* GetName() {
		return "Gaussian Diffusion-Reaction";
	}
};

#endif

