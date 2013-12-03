#ifndef CELLMODEL_H
#define CELLMODEL_H

class CellModel : public BaseIVP {
	FP _cV;  // Cell volume
	FP _bNa; // Bath [Na+]
	FP _bK;  // Bath [K+]
	FP _bCl; // Bath [Cl-]

	FP _R, _T, _F;
	FP _kNa, _kK, _kCl;
	Vec<FP> _P;

	void CalculateCurrents(const FP t, const Vec<FP>& y, FP Ua, FP Ub, Vec<FP>& J) {
		FP expUa = exp(-Ua);
		FP expUb = exp(-Ub);

		// Single ion channels
		J(0) = _P(0) * Ua * (_bNa - y(0)*expUa) / (1-expUa);
		J(1) = _P(1) * Ua * (_bK  - y(1)*expUa) / (1-expUa);
		J(2) = _P(2) * Ua * (_bCl*expUa - y(2)) / (expUa-1);
		J(3) = _P(3) * Ub * (_bK  - y(1)*expUb) / (1-expUb);

		// 3Na/2K ATPase
		FP b = 1 - 6e-3*Ub;
		J(4) = _P(4) * pow(y(0)/(y(0)+_kNa), 3) * pow(_bK/(_bK+_kK), 2) * (6e-3*Ub + b);

		// Na/K/2Cl Cotransporter
		J(5) = _P(5) * (_bNa*_bK*sqr(_bCl) - y(0)*y(1)*sqr(y(2))) / ( (y(0)/_kNa+1) * (y(1)/_kK+1) * sqr(y(2)/_kCl+1) );

		// Na/Cl Cotransporter
		J(6) = _P(6) * (_bNa*_bCl - y(0)*y(2)) / ( (y(0)/_kNa+1) * (y(2)/_kCl+1) );
	}

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		yp.Zero();

		FP Va = -46;
		//FP Vb = -46;

		FP Ua = 1e-3 * Va * _F/(_R*_T);
		//FP Ub = 1e-3 * Vb * _F/(_R*_T);

		// Epsilon for jacobian	
		FP eps = std::numeric_limits<FP>().epsilon();

		Vec<FP> J(7);
		for( short i = 0; ;i++ ) {
			// Begin Jacobian calcs
			FP delta = sqrt(eps*std::max(1e-5, fabs(Ua)));
			CalculateCurrents(t, y, Ua+delta, Ua+delta, J);
			FP jac = (J(0) + J(1) - J(2)) - (J(4) - J(3));
			
			// Current evaluation
			CalculateCurrents(t, y, Ua, Ua, J);
			FP Ia = J(0) + J(1) - J(2);
			FP Ib = J(4) - J(3);
			FP I = Ia - Ib;
			
			// Finish calculating Jacobian
			jac -= I;
			jac /= delta;

			// Newton update
			Ua -= I/jac;
			
			FP absI = fabs(I);
			if( absI < 1e-5 ) {
				//std::cout << "Newton method converged after " << i << " iterations.\n";
				break;
			}

			if( absI > 1e20 )
				throw Exception() << "Newton method blew up.";
			
			if( i > 15 )
				throw Exception() << "Newton method failed to converge.";
		}

		yp(0) = J(0) - 3*J(4) + J(5) + J(6);
		yp(1) = J(1) + J(3) + 2*J(4) + J(5);
		yp(2) = J(2) + 2*J(5) + J(6);
		yp /= _cV;
	}

public:
	CellModel(Hash<ParamValue>& params) : BaseIVP(params) {
		// Cell volume
		_cV = 5e-3;

		// Bath concentrations
		_bNa = 100;
		_bK  = 10;
		_bCl = 100;

		// Channel permeabilities
		_P.Resize(7);
		_P(0) = 2.0e-3; 
		_P(1) = 3.0e-3;
		_P(2) = 1.4e-2;
		_P(3) = 1.2e-2;
		_P(4) = 3.9e+0;
		_P(5) = 1.5e-5;
		_P(6) = 1.0e-4;

		// Saturability values
		_kNa = 3.8;
		_kK  = 7.5;
		_kCl = 26;

		// Thermodynamic constants
		_R = 8.314472;
		_T = 310.15;
		_F = 96485.3365;

		// Initial conditions of the experiment
		_initialCondition.Resize(3);
		_initialCondition[0] = 10; // mmol of Na
		_initialCondition[1] = 80; // mmol of K
		_initialCondition[2] = 10; // mmol of Cl
	}

	IVP_NAME("Epithelial Cell Model")
};

#endif

