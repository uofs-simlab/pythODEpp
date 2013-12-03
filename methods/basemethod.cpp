#include <core/exception.h>
#include <methods/basemethod.h>

BaseMethod::BaseMethod(Hash<ParamValue>& params, BaseIVP* ivp) : _ivp(ivp) {
	_acceptedSteps = 0;
	_dtOld = 0;

	_newtonFail = GetDefaultFP(params, "newton fail", 1e20);
	_newtonTol = GetDefaultFP(params, "newton tol", 1e-8);
	_sparse = (bool)GetDefaultLong(params, "sparse", 0);
	_benchmark = (bool)GetDefaultLong(params, "benchmark", 0);
}

BaseMethod::~BaseMethod() {
}

void BaseMethod::SetSolverVariables(long* acceptedSteps, FP* dtold) {
	_acceptedSteps = acceptedSteps;
	_dtOld = dtold;
}

void BaseMethod::SetAccept(bool accept) {
	_accept = accept;
}

bool BaseMethod::Accept() const {
	return _accept;
}

void BaseMethod::PreStep(const FP tn, FP& dt, Vec<FP>& yn) {
}

void BaseMethod::PostStep(const FP tn, FP& dt, Vec<FP>& yn) {
}

void BaseMethod::UpdateTimestep() {
}

FP BaseMethod::CalcEpsilon(FP tn, FP dt, const Vec<FP>& yn, const Vec<FP>& ynew, FP atol, FP rtol) {
	throw Exception() << "Error calculation is not implemented for " << GetName() << ".";
}

void BaseMethod::GetStats(Hash<ParamValue>& params) const {
}

long BaseMethod::GetAuxOrder() const {
	return GetOrder();
}

