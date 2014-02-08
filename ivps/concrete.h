#ifndef CONCRETE_H
#define CONCRETE_H

#include <ivps/splitivp.h>

#ifdef USE_ADOL_C

Vec<adouble> VecReal(const Vec< std::complex<adouble> >& cv) {
    Vec<adouble> ret(cv.Size());
    for( long i = 0; i < cv.Size(); i++ )
        ret(i) = cv(i).real();
    return ret;
}

#endif

class ConcreteRewetting : public TwoSplittingIVP {
	Vec<FP> _theta;
	Vec<FP> _Ca;
	Vec<FP> _Cb;
	Vec<FP> _Cq;
	Vec<FP> _Cg;

	Vec<FP> _thetaHalf;
	Vec<FP> _D_w_half;
	Vec<FP> _u_half;

	FP _hx;
	FP _rho_csh;
	FP _A, _B;
	FP _phi_pref;
	FP _theta_min, _theta_max;
	FP _theta_rx;
	FP _DCa, _DCb, _DCq;
	FP _k_ads, _k_c3s, _k_c2s, _k_des, _n_c3s, _n_c2s;
	FP _m_csh, _m_c3s, _m_c2s, _m_sol, _rho_sol;
	FP _stoich;
	FP _Xfront0;

	FP _Ca_init;
	FP _Cb_init;
	FP _Cq_init;
	FP _Cg_init;

	long _N;
	long _iMix;
	bool _sinkBC;
	bool _isopropanol;

	void BuildConcreteParameters() {
		FP kfac = 22;

		_k_ads = kfac*1.4640;
		_k_c3s = kfac*1.0109;
		_k_c2s = kfac*0.1382;
		_k_des = 0.0;
		_n_c3s = 2.65;
		_n_c2s = 3.10;
		_stoich = 5.;
	
		// Densities and molar masses
		_rho_csh  = 2.6;
		FP _rho_cem0 = 3.15;
		FP _rho_agg  = 2.60;
		FP _rho_fly  = 2.08;
		FP _rho_sf   = 2.16;
		FP _rho_h2o  = 1.00;
		FP _rho_iso  = 0.785;
	
		// Elemental molar masses [g/mol]:
		FP _m_ca   = 40.078;    // calcium
		FP _m_si   = 28.086;    // silicon
		FP _m_o    = 15.999;    // oxygen
		FP _m_h    =  1.008;    // hydrogen
		FP _m_c    = 12.011;    // carbon
		FP _m_al   = 26.982;    // aluminum
		FP _m_fe   = 55.845;    // iron
	
		// Calculate molar masses for compounds [g/mol]:
		FP _m_h2o   = 2*_m_h + _m_o;                   // water (18.02)
		FP _m_cao   = _m_ca + _m_o;                    // calcium oxide (lime)
		FP _m_sio2  = _m_si + 2*_m_o;                  // silicon dioxide
		FP _m_al2o3 = 2*_m_al + 3*_m_o;                // aluminum oxide
		FP _m_fe2o3 = 2*_m_fe + 3*_m_o;                // ferric oxide
		_m_csh      = 3*_m_cao + 2*_m_sio2 + 3*_m_h2o; // C_3 S_2 H_3 (342.4)
		_m_c3s      = 3*_m_cao + _m_sio2;              // alite,  C_3 S (228.3)
		_m_c2s      = 2*_m_cao + _m_sio2;              // belite, C_2 S (172.2)
		FP _m_c3a   = 3*_m_cao + _m_al2o3;             // tricalcium aluminate, C_3 A (270.2)
		FP _m_c4af  = 4*_m_cao + _m_al2o3 + _m_fe2o3;  // tetracalcium aluminoferrite, C_4 AF (486.0)
		FP _m_iso   = 3*_m_c + 7*_m_h + (_m_o+_m_h);   // isopropanol/C_3 H_7 OH (60.1)

		Vec<FP> wt0[7];	
		for( long i = 0; i < 7; i++ )
			wt0[i].Resize(7);

		wt0[0][0] = 197; wt0[0][1] = 329; wt0[0][2] = 0;   wt0[0][3] = 0;  wt0[0][4] = 921;  wt0[0][5] = 853; wt0[0][6] = 0;
		wt0[1][0] = 182; wt0[1][1] = 250; wt0[1][2] = 250; wt0[1][3] = 0;  wt0[1][4] = 871;  wt0[1][5] = 692; wt0[1][6] = 0;
		wt0[2][0] = 182; wt0[2][1] = 382; wt0[2][2] = 164; wt0[2][3] = 0;  wt0[2][4] = 871;  wt0[2][5] = 692; wt0[2][6] = 0;
		wt0[3][0] = 160; wt0[3][1] = 495; wt0[3][2] = 0;   wt0[3][3] = 43; wt0[3][4] = 1040; wt0[3][5] = 640; wt0[3][6] = 0;
		wt0[4][0] = 197; wt0[4][1] = 329; wt0[4][2] = 0;   wt0[4][3] = 0;  wt0[4][4] = 737;  wt0[4][5] = 853; wt0[4][6] = 82;
		wt0[5][0] = 182; wt0[5][1] = 382; wt0[5][2] = 164; wt0[5][3] = 0;  wt0[5][4] = 697;  wt0[5][5] = 692; wt0[5][6] = 78;
		wt0[6][0] = 160; wt0[6][1] = 495; wt0[6][2] = 0;   wt0[6][3] = 43; wt0[6][4] = 832;  wt0[6][5] = 640; wt0[6][6] = 98;
		
		for( long i = 0; i < 7; i++ )
			wt0[i] *= 0.001;
	
		long _ilast = 6;
		FP _wtcem = wt0[_iMix-1].Splice(1,4).Sum();
		Vec<FP> _wtpc = wt0[_iMix-1].Splice(0,_ilast) / _wtcem;
	
		FP _R_wc      = wt0[_iMix-1][0] / _wtcem;
		FP _R_ac      = wt0[_iMix-1].Splice(4,6).Sum() / _wtcem;
	
		FP _rho_cem   = _rho_cem0 * _wtpc[1] + _rho_fly*_wtpc[2] + _rho_sf*_wtpc[3];
	
		// Hydration fractions
		FP _f_c3s     = 0.60;     // alite  (0.70)
		FP _f_c2s     = 0.20;     // belite (0.25)
		FP _f_c3a     = 0.72;     // C3A
		FP _f_c4af    = 0.59;     // C4AF
	
		// Weight fractions
		FP _omega_c3s = 0.62+0.03; // C3S
		FP _omega_c2s = 0.17;      // C2S
		FP _omega_c3a = 0.11;      // C3A
		FP _omega_c4af = 0.00;     // C4AF
	
		// Derived quantities
		FP _dV_c3s  =  53.28/_m_c3s;  // C3S (hydration)
		FP _dV_c2s  =  39.35/_m_c2s;  // C2S (hydration)
		FP _dV_c3a  = 149.83/_m_c3a;  // C3A
		FP _dV_c4af = 230.00/_m_c4af; // C4AF
	
		if( _isopropanol ) { // (no reactions)
			_rho_sol = _rho_iso;
			_m_sol   = _m_iso;
			_k_c3s   = 0.0;
			_k_c2s   = 0.0;
			_k_ads   = 0.0;
			_k_des   = 0.0;
		} else {
			_rho_sol = _rho_h2o; // solute is water
			_m_sol   = _m_h2o;
		}
	
		// Compute the initial alite and belite concentrations.
		FP _Cmix   = _rho_cem / (_R_wc*_rho_cem/_rho_h2o + _R_ac*_rho_cem/_rho_agg + 1);
		_Ca_init   = (1.0 - _f_c3s) * _omega_c3s * _Cmix; // alite  [g/cm^3]
		_Cb_init   = (1.0 - _f_c2s) * _omega_c2s * _Cmix; // belite [g/cm^3]
		_Cq_init   = 0.0;      // initial CSH(aq) concentration  [g/cm^3]
		_Cg_init   = 0.0;      // initial CSH(gel) concentration [g/cm^3]
	
		_theta_min = 0.04;
		_phi_pref = _Cmix * (_R_wc/_rho_h2o) - _theta_min - _Cmix * \
				(_f_c3s*_omega_c3s*_dV_c3s + _f_c2s*_omega_c2s*_dV_c2s + \
				_f_c3a*_omega_c3a*_dV_c3a + _f_c4af*_omega_c4af*_dV_c4af);
	
		_theta_rx  = _theta_min;
		_theta_max = _phi_pref;
	
		_A         = 0.0028; // water diffusion coeff [cm^2/day]
		_B         = 100;    // water diffusivity exponent
	
		_DCa       = 0.01;   // alite diffusivity [cm^2/day]
		_DCb       = _DCa;   // belite diffusivity
		_DCq       = _DCa;   // CSH(aq) diffusivity
	
		// Grid stuff
		FP _H  = 10.; // Height (in cm)
		_hx = _H/_N;
		_Xfront0 = 0;

		// Dependant variables
		_theta.Resize(_N+2);
		_Ca.Resize(_N+2);
		_Cb.Resize(_N+2);
		_Cq.Resize(_N+2);
		_Cg.Resize(_N+2);

		//InitADOLCParameters();
	}

	void CalculateCommon(const FP t, const Vec<FP>& y) {
		if( _freezeCommon )
			return;
	
		_theta.AssignSplice(y.Splice(0, _N), 1);
	
		_Cg.AssignSplice(y.Splice(4*_N,5*_N) / _theta.Splice(1,_N+1), 1);
		_Cg[0]    = 2*_Cg[1]    - _Cg[2];
		_Cg[_N-1] = 2*_Cg[_N-2] - _Cg[_N-3];
		Vec<FP> CgHalf = (_Cg.Splice(0,_N+1) + _Cg.Splice(1,_N+2)) / 2;

		Vec<FP> phi_p      = _phi_pref * (1.0 - _Cg/(_phi_pref*_rho_csh)).Maximum(0);
		Vec<FP> phi_p_half = _phi_pref * (1 - CgHalf/(_phi_pref*_rho_csh)).Maximum(0);
	
		_theta[0]    = 2*_theta_max - _theta[1] - 2*CgHalf[0]/_rho_csh;
	    _theta[_N+1] = _theta[_N];
		_thetaHalf   = (_theta.Splice(0,_N+1) + _theta.Splice(1,_N+2)) / 2;
	
		// Compute some other cell edge values:
		_D_w_half = _A * VecReal(Vec<CFP>(phi_p_half / _phi_pref).Power(19./6)) * (_B*_thetaHalf).Exp();	

		// Extract remaining solution vectors and set boundary values:
		_Ca.AssignSplice(y.Splice(_N,2*_N) / _theta.Splice(1,_N+1), 1);
		_Ca[_N+1] = _Ca[_N]; // zero Neumann
	
		_Cb.AssignSplice(y.Splice(2*_N,3*_N) / _theta.Splice(1,_N+1), 1);
		_Cb[_N+1] = _Cb[_N]; // zero Neumann
	
		_Cq.AssignSplice(y.Splice(3*_N,4*_N) / _theta.Splice(1,_N+1), 1);
		_Cq[_N+1] = _Cq[_N]; // zero Neumann
	
		_u_half = -_D_w_half * _theta.Diff() / _hx;
		
		// Apply either perfect sink or zero flux BC at x=0.  Perfect sink
		// BC's are implemented only to first order to avoid negative values
		// of concentration.
		if( _sinkBC ) {
			long sinkfac = 0; // 0=first order, -1=second order
			_Ca[0]  = sinkfac*_Ca[1];
			_Cb[0]  = sinkfac*_Cb[1];
			_Cq[0]  = sinkfac*_Cq[1];
		} else {
			FP Aj = _hx * _u_half[0] / (2 * _thetaHalf[0]);
			_Ca[0]   = (1-Aj/_DCa)/(1+Aj/_DCa) * _Ca[1];
			_Cb[0]   = (1-Aj/_DCb)/(1+Aj/_DCb) * _Cb[1];
			_Cq[0]   = (1-Aj/_DCq)/(1+Aj/_DCq) * _Cq[1];
		}
	}

	void Split1(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		CalculateCommon(t, y);
		
		Vec<FP> theta_xx = ( _D_w_half.Splice(1,_N+1)*(_theta.Splice(2,_N+2)-_theta.Splice(1,_N+1)) - 
	                         _D_w_half.Splice(0,_N)*(_theta.Splice(1,_N+1)-_theta.Splice(0,_N)) ) / (_hx*_hx);
		Vec<FP> C_axx = ( _thetaHalf.Splice(1,_N+1) * (_Ca.Splice(2,_N+2)-_Ca.Splice(1,_N+1)) -
	                      _thetaHalf.Splice(0,_N)   * (_Ca.Splice(1,_N+1)-_Ca.Splice(0,_N)) ) / (_hx*_hx);

		Vec<FP> C_bxx = ( _thetaHalf.Splice(1,_N+1) * (_Cb.Splice(2,_N+2)-_Cb.Splice(1,_N+1)) -
	                      _thetaHalf.Splice(0,_N)   * (_Cb.Splice(1,_N+1)-_Cb.Splice(0,_N)) ) / (_hx*_hx);

		Vec<FP> C_qxx = ( _thetaHalf.Splice(1,_N+1) * (_Cq.Splice(2,_N+2)-_Cq.Splice(1,_N+1)) -
	                      _thetaHalf.Splice(0,_N)   * (_Cq.Splice(1,_N+1)-_Cq.Splice(0,_N)) ) / (_hx*_hx);
	
		// Compute the RHS vector
		yp.AssignSplice(theta_xx,              0);
		yp.AssignSplice(_DCa*C_axx,           _N);
		yp.AssignSplice(_DCb*C_bxx,         2*_N);
		yp.AssignSplice(_DCq*C_qxx,         3*_N);
		yp.AssignSplice(Vec<FP>::Zeros(_N), 4*_N);
	}

	void Split2(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		CalculateCommon(t, y);

		// Calculate remaining edge values
		Vec<FP> C_a_half = (_Ca.Splice(0,_N+1) + _Ca.Splice(1,_N+2)) / 2;
		Vec<FP> C_b_half = (_Cb.Splice(0,_N+1) + _Cb.Splice(1,_N+2)) / 2;
		Vec<FP> C_q_half = (_Cq.Splice(0,_N+1) + _Cq.Splice(1,_N+2)) / 2;
	
		Vec<FP> C_ax = ( _u_half.Splice(1,_N+1)*C_a_half.Splice(1,_N+1) - _u_half.Splice(0,_N)*C_a_half.Splice(0,_N) ) / _hx;
		Vec<FP> C_bx = ( _u_half.Splice(1,_N+1)*C_b_half.Splice(1,_N+1) - _u_half.Splice(0,_N)*C_b_half.Splice(0,_N) ) / _hx;
		Vec<FP> C_qx = ( _u_half.Splice(1,_N+1)*C_q_half.Splice(1,_N+1) - _u_half.Splice(0,_N)*C_q_half.Splice(0,_N) ) / _hx;

	    // Reaction terms, written in units of [g/cm^3/day]
		Vec<FP> Rxn_C_a = _k_c3s * VecReal(Vec<CFP>(_Ca.Splice(1,_N+1)).Power(_n_c3s) / pow(CFP(_Ca_init),_n_c3s-1));
		Vec<FP> Rxn_C_b = _k_c2s * VecReal(Vec<CFP>(_Cb.Splice(1,_N+1)).Power(_n_c2s) / pow(CFP(_Cb_init),_n_c2s-1));
		Vec<FP> Rxn_csh = 0.5 * _m_csh * (Rxn_C_a / _m_c3s + Rxn_C_b / _m_c2s);
		
		// This is the saturation factor that multiplies the reaction terms.
		// It should go to zero with the reactive water.
		Vec<FP> thfac = (_theta.Splice(1,_N+1) - _theta_rx).Maximum(0);
		
		// Compute the RHS vector
		yp.AssignSplice(-_stoich*(_m_sol/_m_csh/_rho_sol)*thfac*Rxn_csh,                                  0);
		yp.AssignSplice(-C_ax - thfac*Rxn_C_a,                                                           _N);
		yp.AssignSplice(-C_bx - thfac*Rxn_C_b,                                                         2*_N);
		yp.AssignSplice(-C_qx + thfac*(Rxn_csh-_k_ads*_Cq.Splice(1,_N+1) + _k_des*_Cg.Splice(1,_N+1)), 3*_N);
		yp.AssignSplice(        thfac*(        _k_ads*_Cq.Splice(1,_N+1) - _k_des*_Cg.Splice(1,_N+1)), 4*_N);
	}

public:
	ConcreteRewetting(Hash<ParamValue>& params) : TwoSplittingIVP(params) {
		ParamValue* param;
		if( (param = params.Get("N")) ) _N = param->GetLong();
		else _N = 100;
		if( (param = params.Get("sink bc")) ) _sinkBC = (bool)param->GetLong();
		else _sinkBC = true;
		if( (param = params.Get("isopropanol")) ) _isopropanol = (bool)param->GetLong();
		else _isopropanol = false;
		if( (param = params.Get("iMix")) ) _iMix = param->GetLong();
		else _iMix = 3;

		BuildConcreteParameters();
		Vec<FP> theta0(_N);
		for( long i = 0; i < theta0.Size(); i++ )
			theta0[i] = _hx*(i+0.5) < _Xfront0 ? _theta_max : _theta_min;	

		_initialCondition.Resize(5*_N);
		_initialCondition.AssignSplice(theta0,             0);
		_initialCondition.AssignSplice(_Ca_init*theta0,   _N);
		_initialCondition.AssignSplice(_Cb_init*theta0, 2*_N);
		_initialCondition.AssignSplice(_Cq_init*theta0, 3*_N);
		_initialCondition.AssignSplice(_Cg_init*theta0, 4*_N);
		SetDefaultFP(params,"tf",28.);
	}

	virtual const char* GetName() {
		return "Concrete Rewetting";
	}

// -----------------------------------------------------------------------------
// ADOLC STUFF
//
/*
private:
	Vec<adouble> _thetaADOLC;
	Vec<adouble> _CaADOLC;
	Vec<adouble> _CbADOLC;
	Vec<adouble> _CqADOLC;
	Vec<adouble> _CgADOLC;
	Vec<adouble> _thetaHalfADOLC;
	Vec<adouble> _D_w_halfADOLC;
	Vec<adouble> _u_halfADOLC;

	void InitADOLCParameters() {
		_thetaADOLC.Resize(_N+2);
		_CaADOLC.Resize(_N+2);
		_CbADOLC.Resize(_N+2);
		_CqADOLC.Resize(_N+2);
		_CgADOLC.Resize(_N+2);
	}

	void CalculateCommon(const adouble t, const Vec<adouble>& y) {
		if( _freezeCommon )
			return;
	
		_thetaADOLC.AssignSplice(y.Splice(0, _N), 1);
	
		_CgADOLC.AssignSplice(y.Splice(4*_N,5*_N) / _thetaADOLC.Splice(1,_N+1), 1);
		_CgADOLC[0]    = 2*_CgADOLC[1]    - _CgADOLC[2];
		_CgADOLC[_N-1] = 2*_CgADOLC[_N-2] - _CgADOLC[_N-3];
		Vec<adouble> CgHalf = (_CgADOLC.Splice(0,_N+1) + _CgADOLC.Splice(1,_N+2)) / 2;

		Vec<adouble> phi_p      = _phi_pref * (1.0 - _CgADOLC/(_phi_pref*_rho_csh)).Maximum(0);
		Vec<adouble> phi_p_half = _phi_pref * (1 - CgHalf/(_phi_pref*_rho_csh)).Maximum(0);
	
		_thetaADOLC[0]    = 2*_theta_max - _thetaADOLC[1] - 2*CgHalf[0]/_rho_csh;
	    _thetaADOLC[_N+1] = _thetaADOLC[_N];
		_thetaHalfADOLC   = (_thetaADOLC.Splice(0,_N+1) + _thetaADOLC.Splice(1,_N+2)) / 2;
	
		// Compute some other cell edge values:
		_D_w_halfADOLC = _A * (phi_p_half / _phi_pref).RealPower(19./6) * (_B*_thetaHalfADOLC).Exp();	

		// Extract remaining solution vectors and set boundary values:
		_CaADOLC.AssignSplice(y.Splice(_N,2*_N) / _thetaADOLC.Splice(1,_N+1), 1);
		_CaADOLC[_N+1] = _CaADOLC[_N]; // zero Neumann
	
		_CbADOLC.AssignSplice(y.Splice(2*_N,3*_N) / _thetaADOLC.Splice(1,_N+1), 1);
		_CbADOLC[_N+1] = _CbADOLC[_N]; // zero Neumann
	
		_CqADOLC.AssignSplice(y.Splice(3*_N,4*_N) / _thetaADOLC.Splice(1,_N+1), 1);
		_CqADOLC[_N+1] = _CqADOLC[_N]; // zero Neumann
	
		_u_halfADOLC = -_D_w_halfADOLC * _thetaADOLC.Diff() / _hx;
		
		// Apply either perfect sink or zero flux BC at x=0.  Perfect sink
		// BC's are implemented only to first order to avoid negative values
		// of concentration.
		if( _sinkBC ) {
			long sinkfac = 0; // 0=first order, -1=second order
			_CaADOLC[0]  = sinkfac*_CaADOLC[1];
			_CbADOLC[0]  = sinkfac*_CbADOLC[1];
			_CqADOLC[0]  = sinkfac*_CqADOLC[1];
		} else {
			adouble Aj = _hx * _u_halfADOLC[0] / (2 * _thetaHalfADOLC[0]);
			_CaADOLC[0]   = (1-Aj/_DCa)/(1+Aj/_DCa) * _CaADOLC[1];
			_CbADOLC[0]   = (1-Aj/_DCb)/(1+Aj/_DCb) * _CbADOLC[1];
			_CqADOLC[0]   = (1-Aj/_DCq)/(1+Aj/_DCq) * _CqADOLC[1];
		}
	}

	void Split1(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) {
		CalculateCommon(t, y);
		
		Vec<adouble> theta_xx = ( _D_w_halfADOLC.Splice(1,_N+1)*(_thetaADOLC.Splice(2,_N+2)-_thetaADOLC.Splice(1,_N+1)) - 
	                         _D_w_halfADOLC.Splice(0,_N)*(_thetaADOLC.Splice(1,_N+1)-_thetaADOLC.Splice(0,_N)) ) / (_hx*_hx);
		Vec<adouble> C_axx = ( _thetaHalfADOLC.Splice(1,_N+1) * (_CaADOLC.Splice(2,_N+2)-_CaADOLC.Splice(1,_N+1)) -
	                      _thetaHalfADOLC.Splice(0,_N)   * (_CaADOLC.Splice(1,_N+1)-_CaADOLC.Splice(0,_N)) ) / (_hx*_hx);

		Vec<adouble> C_bxx = ( _thetaHalfADOLC.Splice(1,_N+1) * (_CbADOLC.Splice(2,_N+2)-_CbADOLC.Splice(1,_N+1)) -
	                      _thetaHalfADOLC.Splice(0,_N)   * (_CbADOLC.Splice(1,_N+1)-_CbADOLC.Splice(0,_N)) ) / (_hx*_hx);

		Vec<adouble> C_qxx = ( _thetaHalfADOLC.Splice(1,_N+1) * (_CqADOLC.Splice(2,_N+2)-_CqADOLC.Splice(1,_N+1)) -
	                      _thetaHalfADOLC.Splice(0,_N)   * (_CqADOLC.Splice(1,_N+1)-_CqADOLC.Splice(0,_N)) ) / (_hx*_hx);
	
		// Compute the RHS vector
		yp.AssignSplice(theta_xx,              0);
		yp.AssignSplice(_DCa*C_axx,           _N);
		yp.AssignSplice(_DCb*C_bxx,         2*_N);
		yp.AssignSplice(_DCq*C_qxx,         3*_N);
		yp.AssignSplice(Vec<adouble>::Zeros(_N), 4*_N);
	}

	void Split2(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) {
		CalculateCommon(t, y);

		// Calculate remaining edge values
		Vec<adouble> C_a_half = (_CaADOLC.Splice(0,_N+1) + _CaADOLC.Splice(1,_N+2)) / 2;
		Vec<adouble> C_b_half = (_CbADOLC.Splice(0,_N+1) + _CbADOLC.Splice(1,_N+2)) / 2;
		Vec<adouble> C_q_half = (_CqADOLC.Splice(0,_N+1) + _CqADOLC.Splice(1,_N+2)) / 2;
	
		Vec<adouble> C_ax = ( _u_halfADOLC.Splice(1,_N+1)*C_a_half.Splice(1,_N+1) - _u_halfADOLC.Splice(0,_N)*C_a_half.Splice(0,_N) ) / _hx;
		Vec<adouble> C_bx = ( _u_halfADOLC.Splice(1,_N+1)*C_b_half.Splice(1,_N+1) - _u_halfADOLC.Splice(0,_N)*C_b_half.Splice(0,_N) ) / _hx;
		Vec<adouble> C_qx = ( _u_halfADOLC.Splice(1,_N+1)*C_q_half.Splice(1,_N+1) - _u_halfADOLC.Splice(0,_N)*C_q_half.Splice(0,_N) ) / _hx;

	    // Reaction terms, written in units of [g/cm^3/day]
		Vec<adouble> Rxn_C_a;// = _k_c3s * VecReal(Vec< std::complex<adouble> >(_CaADOLC.Splice(1,_N+1)).Power(adouble(_n_c3s)) / pow(std::complex<adouble>(_Ca_init),_n_c3s-1));
		Vec<adouble> Rxn_C_b;// = _k_c2s * VecReal(Vec< std::complex<adouble> >(_CbADOLC.Splice(1,_N+1)).Power(adouble(_n_c2s)) / pow(std::complex<adouble>(_Cb_init),_n_c2s-1));
		Vec<adouble> Rxn_csh = 0.5 * _m_csh * (Rxn_C_a / _m_c3s + Rxn_C_b / _m_c2s);
		
		// This is the saturation factor that multiplies the reaction terms.
		// It should go to zero with the reactive water.
		Vec<adouble> thfac = (_thetaADOLC.Splice(1,_N+1) - _theta_rx).Maximum(0);
		
		// Compute the RHS vector
		yp.AssignSplice(-_stoich*(_m_sol/_m_csh/_rho_sol)*thfac*Rxn_csh,                                  0);
		yp.AssignSplice(-C_ax - thfac*Rxn_C_a,                                                           _N);
		yp.AssignSplice(-C_bx - thfac*Rxn_C_b,                                                         2*_N);
		yp.AssignSplice(-C_qx + thfac*(Rxn_csh-_k_ads*_CqADOLC.Splice(1,_N+1) + _k_des*_CgADOLC.Splice(1,_N+1)), 3*_N);
		yp.AssignSplice(        thfac*(        _k_ads*_CqADOLC.Splice(1,_N+1) - _k_des*_CgADOLC.Splice(1,_N+1)), 4*_N);
	}*/
};

#endif

