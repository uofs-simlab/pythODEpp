#ifndef NONSTIFF_C5_H
#define NONSTIFF_C5_H

class NonstiffC5 : public BaseIVP {
protected:
	FP _k2;
	FP _m0;
	FP _m[5];

	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
		FP ypos[3][5];
		FP ypp[3][5];
		FP rSquared[5];
		FP rCubed[5];
		FP dSquared[5][5];
		FP dCubed[5][5];

		for( short j = 0; j < 5; j++ )
			for( short i = 0; i < 3; i++ )
				ypos[i][j] = y[i+j*3];

		for( short j = 0; j < 5; j++ ) {
			rSquared[j] = sqr(ypos[0][j]) + sqr(ypos[1][j]) + sqr(ypos[2][j]);
			rCubed[j] = rSquared[j]*sqrt(rSquared[j]);
		}

		for( short k = 0; k < 5; k++ ) {
			for( short j = 0; j < 5; j++ ) {
				dSquared[k][j] = sqr(ypos[0][k]-ypos[0][j]) + sqr(ypos[1][k]-ypos[1][j]) + sqr(ypos[2][k]-ypos[2][j]);
				dCubed[k][j] = dSquared[k][j]*sqrt(dSquared[k][j]);
			}
		}

		for( short i = 0; i < 3; i++ ) {
			for( short j = 0; j < 5; j++ ) {
				ypp[i][j] = -(_m0 + _m[j]) * ypos[i][j]/rCubed[j];
				for( short k = 0; k < 5; k++ ) {
					if( k != j )
						ypp[i][j] += _m[k]*(((ypos[i][k]-ypos[i][j])/dCubed[j][k]) - ypos[i][k]/rCubed[k]);
				}
				ypp[i][j] *= _k2;
			}
		}

		for( short i = 0; i < 15; i++ )
			yp[i] = y[15+i];
		for( short j = 0; j < 3; j++ )
			for( short k = 0; k < 5; k++ )
				yp[15+j+k*3] = ypp[j][k];
	}

public:
	NonstiffC5(Hash<ParamValue>& params) : BaseIVP(params) {
		params["tf"].SetFP(20);
		
		_k2 = 2.9591220828;
		_m0 = 1.0000059768;
		_m[0] = 0.00095478610404;
		_m[1] = 0.00028558373315;
		_m[2] = 0.000043727316454;
		_m[3] = 0.000051775913844;
		_m[4] = 0.000002777777777;

		_initialCondition.Resize(30);
		_initialCondition[0]  = 3.42947415189e0;
		_initialCondition[1]  = 3.35386959711e0;
		_initialCondition[2]  = 1.35494901715e0;
		_initialCondition[3]  = 6.64145542550e0;
		_initialCondition[4]  = 5.97156957878e0;
		_initialCondition[5]  = 2.18231499728e0;
		_initialCondition[6]  = 11.2630437207e0;
		_initialCondition[7]  = 14.6952576794e0;
		_initialCondition[8]  = 6.27960525067e0;
		_initialCondition[9]  = -30.1552268759e0;
		_initialCondition[10]  = 1.65699966404e0;
		_initialCondition[11] = 1.43785752721e0;
		_initialCondition[12] = -21.1238353380e0;
		_initialCondition[13] = 28.4465098142e0;
		_initialCondition[14] = 15.3882659679e0;
		_initialCondition[15] = -.557160570446e0;
		_initialCondition[16] = 0.505696783289e0;
		_initialCondition[17] = 0.230578543901e0;
		_initialCondition[18] = -0.415570776342e0;
		_initialCondition[19] = 0.365682722812e0;
		_initialCondition[20] = 0.169143213293e0;
		_initialCondition[21] = -0.325325669158e0;
		_initialCondition[22] = 0.189706021964e0;
		_initialCondition[23] = 0.0877265322780e0;
		_initialCondition[24] = -0.0240476254170e0;
		_initialCondition[25] = -0.287659532608e0;
		_initialCondition[26] = -0.117219543175e0;
		_initialCondition[27] = -0.176860753121e0;
		_initialCondition[28] = -0.216393453025e0;
		_initialCondition[29] = -0.0148647893090e0;
	}
	
	virtual const char* GetName() {
		return "Nonstiff C5";
	}
};

#endif

