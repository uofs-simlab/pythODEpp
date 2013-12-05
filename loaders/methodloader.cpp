#include <methods/erk.h>
#include <methods/irk.h>
#include <methods/ark.h>
#include <methods/rkc.h>
#include <methods/exprk.h>

#define METHODCASE(methodclass) if( method == #methodclass ) return new methodclass(params,ivp);

BaseMethod* AllocMethod(Hash<ParamValue>& params, BaseIVP* ivp) {
	std::string method = params["method"].GetString();

	METHODCASE(ForwardEuler)
	METHODCASE(BackwardEuler)
	METHODCASE(RK4)
	METHODCASE(RK38)
	METHODCASE(Runge2)
	METHODCASE(Runge3)
	METHODCASE(Kutta3)
	METHODCASE(Heun2)
	METHODCASE(Heun3)
	METHODCASE(RKF45)
	METHODCASE(FEHL78)
	METHODCASE(DOPR54)
	METHODCASE(Merson43)
	METHODCASE(Zonneveld43)
	METHODCASE(Verner65)
	METHODCASE(RKC1)
	METHODCASE(RKC2)
	METHODCASE(PRKC)
	METHODCASE(IRKC)
	METHODCASE(ARK1)
	METHODCASE(ARK3)
	METHODCASE(ARK4)
	METHODCASE(ARK5)
	METHODCASE(Radau5)
	METHODCASE(RODAS)
	METHODCASE(DIRKCF1)
	METHODCASE(DIRKCF2)
	METHODCASE(DIRKCF3)
	METHODCASE(ERKCF2)

	throw Exception() << "Method " << method << " has not been defined.";
}

