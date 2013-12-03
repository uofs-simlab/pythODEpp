#include <core/exception.h>
#include <core/timer.h>
#include <ivps/splitivp.h>

void TwoSplittingIVP::FreezeCommon(bool fc) {
	_freezeCommon = fc;
}

FP TwoSplittingIVP::CustomMax0(FP a) {
	return max(a,0.);
}

#ifdef USE_ADOL_C
adouble TwoSplittingIVP::CustomMax0(adouble a) {
	adouble ret = 0;
	condassign(ret, a, a);
	return ret;
}
#endif

FP TwoSplittingIVP::ProcessConditional(FP c, FP t, FP f) {
	return c > 0 ? t : f;
}   

#ifdef USE_ADOL_C
adouble TwoSplittingIVP::ProcessConditional(adouble c, adouble t, adouble f) {
	adouble ret;
	condassign(ret, c, t, f);
	return ret;
}   
#endif

void TwoSplittingIVP::RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) {
	CalculateCommon(t,y);
	FreezeCommon(true);

	Split1(t, y, yp);;

	Vec<FP> split2(y.Size());
	Split2(t, y, split2);
	yp += split2;

	FreezeCommon(false);

	_fEvals++;
	_gEvals++;
}   

void TwoSplittingIVP::PhysicalSplit(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& yp) {
	if( split == 1 ) { 
		Split1(t, y, yp);
		_fEvals++;
	} else if( split == 2 ) { 
		Split2(t,y,yp);
		_gEvals++;
	} else {
		throw Exception() << GetName() << " is two-splitting.";
	}
}

void TwoSplittingIVP::PhysicalSplitMat(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& mat) {
	if( split == 1 ) { 
		SplitMat1(t, y, mat);
	} else if( split == 2 ) { 
		SplitMat2(t, y, mat);
	} else {
		throw Exception() << GetName() << " is two-splitting.";
	}
}

#ifdef USE_ADOL_C
void TwoSplittingIVP::RHS(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) {
	CalculateCommon(t,y);
	FreezeCommon(true);

	Split1(t, y, yp);

	Vec<adouble> split2(y.Size());
	Split2(t, y, split2);
	yp += split2;

	FreezeCommon(false);

	_fEvals++;
	_gEvals++;
}

void TwoSplittingIVP::PhysicalSplit(unsigned short split, const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) {
	if( split == 1 ) { 
		Split1(t, y, yp);
		_fEvals++;
	} else if( split == 2 ) { 
		Split2(t,y,yp);
		_gEvals++;
	} else {
		throw Exception() << GetName() << " is two-splitting.";
	}
}
#endif

void TwoSplittingIVP::SplitMat1(const FP t, const Vec<FP>& y, Mat<FP>& mat) {
	throw Exception() << GetName() << " does not implement a matrix for the first split component.";
}

void TwoSplittingIVP::SplitMat2(const FP t, const Vec<FP>& y, Mat<FP>& mat) {
	throw Exception() << GetName() << " does not implement a matrix for the second split component.";
}

#ifdef USE_ADOL_C
void TwoSplittingIVP::Split1(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) {
	throw Exception() << GetName() << " does not implement a split term 1 for ADOL-C.";
}

void TwoSplittingIVP::Split2(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) {
	throw Exception() << GetName() << " does not implement a split term 2 for ADOL-C.";
}
#endif


#ifdef USE_SUITESPARSE

void TwoSplittingIVP::PhysicalSplitMatSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& mat) {
	if( split == 1 ) { 
		SplitMat1Sparse(t, y, mat);
	} else if( split == 2 ) { 
		SplitMat2Sparse(t, y, mat);
	} else {
		throw Exception() << GetName() << " is two-splitting.";
	}
}

void TwoSplittingIVP::SplitMat1Sparse(const FP t, const Vec<FP>& y, CSRMat<FP>& mat) {
	throw Exception() << GetName() << " does not implement a sparse matrix for the first split component.";
}

void TwoSplittingIVP::SplitMat2Sparse(const FP t, const Vec<FP>& y, CSRMat<FP>& mat) {
	throw Exception() << GetName() << " does not implement a sparse matrix for the second split component.";
}

#endif

TwoSplittingIVP::TwoSplittingIVP(Hash<ParamValue>& params) : BaseIVP(params,2), _fEvals(0), _gEvals(0), _freezeCommon(false) {
}

void TwoSplittingIVP::GetStats(Hash<ParamValue>& params) const {
	BaseIVP::GetStats(params);
	params["f evaluations"] = _fEvals;
	params["g evaluations"] = _gEvals;
}

void TwoSplittingIVP::PrintStats() const {
	BaseIVP::PrintStats();
	std::cout << "f evaluations = " << _fEvals << "\n";
	std::cout << "g evaluations = " << _gEvals << "\n";
}

