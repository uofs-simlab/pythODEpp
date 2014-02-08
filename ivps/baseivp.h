// Copyright (c) 2013, Adam Preuss
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the author nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef BASE_IVP_H
#define BASE_IVP_H

#include <core/common.h>
#include <core/hash.h>
#include <core/paramvalue.h>
#include <core/vec.h>
#include <core/mat.h>
#include <core/csrmat.h>
#include <core/timer.h>

#ifdef USE_ADOL_C
	#include <adolc/adolc.h>
	#include <adolc/adolc_sparse.h>
#endif

#define IVP_NAME(name) virtual const char* GetName() { return name; } 
#define SPLIT_FP(fname,target) void fname(const FP t, const Vec<FP>& y, Vec<FP>& yp) { target(t,y,yp); }

#ifdef USE_ADOL_C
	#define SPLIT_ADOLC(fname,target) void fname(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) { target(t,y,yp); }
#else
	#define SPLIT_ADOLC(fname,target)
#endif

class BaseIVP {
protected:
	// Approaches to calculating derivatives and Jacobians
	enum DType {
		D_ANALYTIC = 0,
		D_AUTODIFF = 1,
		D_FORWARD = 2,
		D_CENTRED = 3
	};

	// Standard IVP parameters
	FP _initialTime;
	FP _finalTime;
	Vec<FP> _initialCondition;

	// Derivative and Jacobian parameters
	FP _dtDelta;
	DType _dtType;
	FP _jacDelta;
	DType _jacType;
	bool _jacFrozen;
	bool _jacSplitting;
	FP _jacScaling;

	// Storage locations of the Jacobian matrices
private:
	unsigned short _splitCount;
protected:
	BaseMat<FP>** _splitMats;
	BaseMat<FP>** _splitJacs;

	// Finite difference order
	long _fdorder;

	virtual void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac);

	void JacAutodiff(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac);
	void JacForward(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac);
	void JacCentred(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac);	

	void DtAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& dfdt);
	void DtAutodiff(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& dfdt);
	void DtForward(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& dfdt);
	void DtCentred(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& dfdt);

	virtual void PhysicalSplit(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& yp);
	virtual void PhysicalSplitMat(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& mat);
	virtual void RHS(const FP t, const Vec<FP>& y, Vec<FP>& yp) = 0;
	
	virtual void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac);

	void JacAutodiffSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac);
	void JacForwardSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac);
	void JacCentredSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac);

	virtual void PhysicalSplitMatSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& mat);

#ifdef USE_ADOL_C
	virtual void RHS(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp);
	virtual void PhysicalSplit(unsigned short split, const adouble t, const Vec<adouble>& y, Vec<adouble>& yp);
#endif

public:
	BaseIVP(Hash<ParamValue>& params, unsigned short splitting = 0);
	virtual ~BaseIVP();
	
	void InitializeDerivatives();

	void GetInitialCondition(Vec<FP>& ic);

	void FreezeJacobian(bool jf);
	bool JacobianSplitting() const;

	void operator()(const FP t, const Vec<FP>& y, Vec<FP>& yp, unsigned short split = 0);
	void RHSTimeDt(const FP t, const Vec<FP>& y, Vec<FP>& pfpt, unsigned short split = 0);

	const BaseMat<FP>* SplitMat(const FP t, const Vec<FP>& y, unsigned short split);
	const BaseMat<FP>* Jac(const FP t, const Vec<FP>& y, unsigned short split = 0);
	const BaseMat<FP>* SplitMatSparse(const FP t, const Vec<FP>& y, unsigned short split);
	const BaseMat<FP>* JacSparse(const FP t, const Vec<FP>& y, unsigned short split = 0);

	virtual void GetStats(Hash<ParamValue>& params) const;
	virtual void PrintStats() const;
	
	long Size() const;
	virtual const char* GetName() = 0;

#ifdef USE_ADOL_C
	static void ToADOLC(const Vec<FP>& vec, Vec<adouble>& adolc);
	static void FromADOLC(Vec<adouble>& adolc, Vec<FP>& vec);
#endif
};

#endif

