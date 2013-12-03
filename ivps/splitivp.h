#ifndef SPLIT_IVP_H
#define SPLIT_IVP_H

#include <ivps/baseivp.h>

#define LINK_TWOSPLIT SPLIT_FP(Split1, Split1Internal) \
					  SPLIT_FP(Split2, Split2Internal) \
					  SPLIT_ADOLC(Split1, Split1Internal) \
					  SPLIT_ADOLC(Split2, Split2Internal)

class TwoSplittingIVP : public BaseIVP {
protected:
	long _fEvals;
	long _gEvals;
	bool _freezeCommon;

	void FreezeCommon(bool fc);

	FP CustomMax0(FP a);
	FP ProcessConditional(FP c, FP t, FP f);

#ifdef USE_ADOL_C
	adouble CustomMax0(adouble a);
	adouble ProcessConditional(adouble c, adouble t, adouble f);
#endif
	
	// TODO FIXME Coalesce into a single function
	void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp);
	void PhysicalSplit(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& yp);
	void PhysicalSplitMat(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& mat);
#ifdef USE_ADOL_C
	void RHS(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp);
	void PhysicalSplit(unsigned short split, const adouble t, const Vec<adouble>& y, Vec<adouble>& yp);
#endif

	virtual void SplitMat1(const FP t, const Vec<FP>& y, Mat<FP>& mat);
	virtual void SplitMat2(const FP t, const Vec<FP>& y, Mat<FP>& mat);
	
	virtual void Split1(const FP t, const Vec<FP>& y, Vec<FP>& yp) = 0;
	virtual void Split2(const FP t, const Vec<FP>& y, Vec<FP>& yp) = 0;
	virtual void CalculateCommon(const FP t, const Vec<FP>& y) { }
#ifdef USE_ADOL_C
	virtual void Split1(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp);
	virtual void Split2(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp);
	virtual void CalculateCommon(const adouble t, const Vec<adouble>& y) { }
#endif

#ifdef USE_SUITESPARSE
	void PhysicalSplitMatSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& mat);
	virtual void SplitMat1Sparse(const FP t, const Vec<FP>& y, CSRMat<FP>& mat);
	virtual void SplitMat2Sparse(const FP t, const Vec<FP>& y, CSRMat<FP>& mat);	
#endif

public:
	TwoSplittingIVP(Hash<ParamValue>& params);

	virtual void GetStats(Hash<ParamValue>& params) const;
	virtual void PrintStats() const;
};

#endif

