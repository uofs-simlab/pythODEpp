#include <ivps/baseivp.h>
#include <ivps/bari.h>
#include <ivps/zbinden.h>
#include <ivps/nonstiff.h>
#include <ivps/concrete.h>
#include <ivps/cellmodel.h>
#include <ivps/plate.h>
#include <ivps/cusp.h>
#include <ivps/brusselator2d.h>
#include <ivps/burgers.h>
#include <ivps/chemotaxis.h>
#include <ivps/rkp.h>
#include <ivps/combustionard.h>
#include <ivps/scottwangshowalter.h>
#include <ivps/gaussian.h>
#include <ivps/heattransfer.h>
#include <ivps/allencahn.h>

#include <ivps/AllenCahn2DF.h>
#include <ivps/VanDerPol.h>
#include <ivps/ViscousBurgers.h>

#define IVPCASE(ivpclass) if( ivp == #ivpclass ) return new ivpclass(params);

BaseIVP* AllocIVP(Hash<ParamValue>& params) {
	std::string ivp = params["ivp"].GetString();

	IVPCASE(Beam)
	IVPCASE(VDPOL)
	IVPCASE(ConcreteRewetting)
	IVPCASE(CellModel)
	IVPCASE(PLATE)
	IVPCASE(CUSP)
	IVPCASE(AdvectionDiffusion1D)
	IVPCASE(Brusselator1D)
	IVPCASE(Brusselator2D)
	IVPCASE(PIDE)
	IVPCASE(Burgers)
	IVPCASE(Angiogenesis1D)
	IVPCASE(RKP)
	IVPCASE(CombustionARD)
	IVPCASE(ScottWangShowalter)
	IVPCASE(Gaussian)
	IVPCASE(HeatTransfer)
	IVPCASE(AllenCahn)

	IVPCASE(NonstiffA1) IVPCASE(NonstiffA2) IVPCASE(NonstiffA3) IVPCASE(NonstiffA4) IVPCASE(NonstiffA5)
	IVPCASE(NonstiffB1) IVPCASE(NonstiffB2) IVPCASE(NonstiffB3) IVPCASE(NonstiffB4) IVPCASE(NonstiffB5)
	IVPCASE(NonstiffC1) IVPCASE(NonstiffC2) IVPCASE(NonstiffC3) IVPCASE(NonstiffC4) IVPCASE(NonstiffC5)
	IVPCASE(NonstiffD1) IVPCASE(NonstiffD2) IVPCASE(NonstiffD3) IVPCASE(NonstiffD4) IVPCASE(NonstiffD5)
	IVPCASE(NonstiffE1) IVPCASE(NonstiffE2) IVPCASE(NonstiffE3) IVPCASE(NonstiffE4) IVPCASE(NonstiffE5)
	IVPCASE(NonstiffF1) IVPCASE(NonstiffF2) IVPCASE(NonstiffF3) IVPCASE(NonstiffF4) IVPCASE(NonstiffF5)

	IVPCASE(AllenCahn2DF)
	IVPCASE(VanDerPol)
	IVPCASE(ViscousBurgers)
    
	throw Exception() << "IVP " << ivp << " has not been defined.";
}

