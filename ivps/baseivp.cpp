#include <core/exception.h>
#include <ivps/baseivp.h>

// -----------------------------------------------------------------------------
// Code for calculating Jacobians and time derivatives.
//
void BaseIVP::JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
	throw Exception() << "Analytic Jacobian is not defined for " << GetName() << ".";
}

void BaseIVP::JacAutodiff(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
#ifdef USE_ADOL_C
	long n = y.Size();

	Vec<FP> in(n+1);
	in.AssignSplice(y, 0);
	in[n] = t;

	double** J = myalloc2(n+1,n);
	jacobian(split, n, n+1, *in, J);
	for( long i = 0; i < n; i++ )
		for( long j = 0; j < n; j++ )
			jac(i,j) = J[i][j];
	myfree2(J);
#else
	throw Exception() << "AD Jacobian is unavailable for " << GetName() << ". Recompile with ADOL-C support.";
#endif
}

void BaseIVP::JacForward(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
	FP eps = std::numeric_limits<FP>().epsilon();

	Vec<FP> f1(y.Size()), f2(y.Size());
	(*this)(t, y, f1, split);
	
	Vec<FP> offset = y;
	for( long j = 0; j < y.Size(); j++ ) {
		FP delta = sqrt(eps*std::max(_jacDelta, fabs(y[j])));
		// Perturb index
		offset(j) += delta;
		(*this)(t, offset, f2, split);
		for( long i = 0; i < y.Size(); i++ )
			jac(i,j) = (f2(i)-f1(i)) / delta;
		// Restore twiddled index
		offset(j) = y(j);
	}
}

void BaseIVP::JacCentred(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac) {
	FP eps = std::numeric_limits<FP>().epsilon();

	Vec<FP> f1(y.Size()), f2(y.Size());	
	Vec<FP> offset1 = y;
	Vec<FP> offset2 = y;
	for( long j = 0; j < y.Size(); j++ ) {
		FP delta = sqrt(eps*std::max(_jacDelta, fabs(y[j])));
		// Perturb indices
		offset1(j) -= delta;
		offset2(j) += delta;
		(*this)(t, offset1, f1, split);
		(*this)(t, offset2, f2, split);
		for( long i = 0; i < y.Size(); i++ )
			jac(i,j) = (f2(i)-f1(i)) / (2*delta);
		// Restore twiddled indices
		offset1(j) = y(j);
		offset2(j) = y(j);
	}
}

void BaseIVP::DtAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& pfpt) {
	throw Exception() << "Analytic time derivative is not defined for " << GetName() << ".";
}

void BaseIVP::DtAutodiff(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& pfpt) {
#ifdef USE_ADOL_C
	throw Exception() << "Analytic time derivative is not yet available in ADOL-C.";
#else
	throw Exception() << "Analytic time derivative is unavailable for " << GetName() << ". Recompile with ADOL-C support.";
#endif
}

void BaseIVP::DtForward(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& pfpt) {
	FP eps = std::numeric_limits<FP>().epsilon();
	FP delta = sqrt(eps*std::max(_dtDelta, fabs(t)));

	Vec<FP> f1(y.Size());
	(*this)(t, y, f1, split);
	(*this)(t+delta, y, pfpt, split);
	pfpt -= f1;
	pfpt /= delta;
}

void BaseIVP::DtCentred(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& pfpt) {
	FP eps = std::numeric_limits<FP>().epsilon();
	FP delta = sqrt(eps*std::max(_dtDelta, fabs(t)));

	Vec<FP> f1(y.Size());
	(*this)(t-delta, y, f1, split);
	(*this)(t+delta, y, pfpt, split);
	pfpt -= f1;
	pfpt /= (2*delta);
}

// -----------------------------------------------------------------------------
// Templates for splittings.
//
void BaseIVP::PhysicalSplit(unsigned short split, const FP t, const Vec<FP>& y, Vec<FP>& yp) {
	throw Exception() << "Physical splitting is not implemented for " << GetName() << ".";
}

void BaseIVP::PhysicalSplitMat(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& mat) {
	throw Exception() << GetName() << " is either not split or does not provide matrices for its splitting.";
}

// -----------------------------------------------------------------------------
// Definitions for sparsity, which are almost identical to the above definitions
//
void BaseIVP::JacForwardSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
	Mat<FP> m(y.Size(), y.Size());
	JacForward(split, t, y, m);
	jac = m;
}

void BaseIVP::JacCentredSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
	Mat<FP> m(y.Size(), y.Size());
	JacCentred(split, t, y, m);
	jac = m;
}

void BaseIVP::JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
	throw Exception() << "Sparse analytic Jacobian is not defined for " << GetName() << ".";
}

void BaseIVP::JacAutodiffSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& jac) {
#ifdef USE_ADOL_C
	long n = y.Size();

	Vec<FP> in(n+1);
	in.AssignSplice(y, 0);
	in[n] = t;

	unsigned int* rind = 0;
	unsigned int* cind = 0;
	double* values = 0;
	int count;
	int options[] = { 0, 0, 0, 0 };

	sparse_jac(split, n, n+1, 0, *in, &count, &rind, &cind, &values, options);

	FP* csrData = new FP[count];
	long* csrCols = new long[count];
	long* csrRows = new long[n];

	long realCount = count;
	for( long i = 0, j = 0, r = 0; i < count; i++ ) {
		if( i == 0 || rind[i] != rind[i-1] ) {
			csrRows[r++] = j;
		}

		if( cind[i] >= n ) {
			realCount--;
			continue;
		}

		csrData[j] = values[i];
		csrCols[j] = cind[i];
		j++;
	}

	jac = CSRMat<FP>(csrData, csrCols, csrRows, n, n, realCount);

	delete [] csrData;
	delete [] csrCols;
	delete [] csrRows;

	free(rind);
	free(cind);
	free(values);
#else
	throw Exception() << "Sparse AD Jacobian is unavailable for " << GetName() << ". Recompile with ADOL-C support.";
#endif
}

void BaseIVP::PhysicalSplitMatSparse(unsigned short split, const FP t, const Vec<FP>& y, CSRMat<FP>& mat) {
	throw Exception() << GetName() << " is either not split or does not provide sparse matrices for its splitting.";
}

#ifdef USE_ADOL_C
void BaseIVP::RHS(const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) {
	throw Exception() << GetName() << " does not implement RHS for ADOL-C.";
}

void BaseIVP::PhysicalSplit(unsigned short split, const adouble t, const Vec<adouble>& y, Vec<adouble>& yp) {
	throw Exception() << GetName() << " does not implement physical split for ADOL-C.";
}
#endif

BaseIVP::BaseIVP(Hash<ParamValue>& params, unsigned short splitting) : _initialTime(0), _finalTime(0), _dtDelta(1e-5), _dtType(D_FORWARD), _jacDelta(1e-5), _jacType(D_FORWARD), _jacFrozen(false), _jacSplitting(false), _jacScaling(1.), _splitCount(splitting), _fdorder(2) {
	ParamValue* pv;
	if( (pv = params.Get("jacobian splitting")) )
		_jacSplitting = (bool)pv->GetLong();

	if( (pv = params.Get("jacobian scaling")) && _jacSplitting )
		_jacScaling = pv->GetFP();

	if( (pv = params.Get("jacobian")) ) {
		if( std::string(pv->GetString()) == "Analytic" )
			_jacType = D_ANALYTIC;
		else if( std::string(pv->GetString()) == "Autodiff" )
			_jacType = D_AUTODIFF;
		else if( std::string(pv->GetString()) == "Forward" )
			_jacType = D_FORWARD;
		else if( std::string(pv->GetString()) == "Centred" )
			_jacType = D_CENTRED;
		else
			throw Exception() << "Unknown Jacobian type " << pv->GetString() << ".";
	}

	if( (pv = params.Get("dfdt")) ) {
		if( std::string(pv->GetString()) == "Analytic" )
			_dtType = D_ANALYTIC;
		else if( std::string(pv->GetString()) == "Autodiff" )
			_dtType = D_AUTODIFF;
		else if( std::string(pv->GetString()) == "Forward" )
			_dtType = D_FORWARD;
		else if( std::string(pv->GetString()) == "Centred" )
			_dtType = D_CENTRED;
		else
			throw Exception() << "Unknown time derivative type " << pv->GetString() << ".";
	}

	if( (pv = params.Get("dt delta")) )
		_dtDelta = pv->GetFP();
	
	if( (pv = params.Get("jac delta")) )
		_jacDelta = pv->GetFP();
	
	if( (pv = params.Get("fdorder")) )
		_fdorder = pv->GetLong();

	if( _jacSplitting )
		_splitCount = 0;

	_splitMats = new BaseMat<FP>*[_splitCount+1];
	_splitJacs = new BaseMat<FP>*[_splitCount+1];

	for( unsigned short i = 0; i <= _splitCount; i++ ) {
		_splitMats[i] = 0;
		_splitJacs[i] = 0;
	}
}

BaseIVP::~BaseIVP() {
	for( unsigned short i = 0; i <= _splitCount; i++ ) {
		if( _splitMats[i] )
			delete _splitMats[i];
		if( _splitJacs[i] )
			delete _splitMats[i];
	}

	if( _splitMats ) delete [] _splitMats;
	if( _splitJacs ) delete [] _splitJacs;
}

void BaseIVP::InitializeDerivatives() {
	if( _jacType != D_AUTODIFF )
		return;

#ifdef USE_ADOL_C
	Vec<FP> dummyVec(_initialCondition);

	adouble t = 0;
	Vec<adouble> y(_initialCondition);
	Vec<adouble> yp(_initialCondition);

	for( long i = 0; i <= _splitCount; i++ ) {
		trace_on(i);
		ToADOLC(dummyVec, y); 
		t <<= 0;

		if( i == 0 )
			RHS(t, y, yp);
		else
			PhysicalSplit(i, t, y, yp);

		FromADOLC(yp, dummyVec);    
		trace_off();
	}
#endif
}

void BaseIVP::GetInitialCondition(Vec<FP>& ic) {
	ic = _initialCondition;
}

void BaseIVP::FreezeJacobian(bool jf) {
	_jacFrozen = jf;
}

bool BaseIVP::JacobianSplitting() const {
	return _jacSplitting;
}

void BaseIVP::operator()(const FP t, const Vec<FP>& y, Vec<FP>& yp, unsigned short split) {
	// No splitting is an easy case
	if( split == 0 ) {
		RHS(t, y, yp);
		return;
	}

	// Physical splitting is also pretty easy. We can ship it off right away.
	if( !_jacSplitting ) {
		PhysicalSplit(split, t, y, yp);
		return;
	}

	// We are using Jacobian splitting (so presumably the jacobian must be frozen at this point)
	if( !_splitJacs[0] )
		throw Exception() << "Jacobian required by Jacobian splitting has not yet been initialized.";

	if( split == 1 ) {
		_splitJacs[0]->VectorMult(y, yp);
		return;
	}

	if( split == 2 ) {
		(*this)(t,y,yp);
		Vec<FP> temp(y.Size());
		_splitJacs[0]->VectorMult(y,temp);
		yp -= temp;
		return;
	}
	
	throw Exception() << "Jacobian splitting is two-splitting";
}

void BaseIVP::RHSTimeDt(const FP t, const Vec<FP>& y, Vec<FP>& pfpt, unsigned short split) {
	switch( _dtType ) {
	case D_ANALYTIC:
		DtAnalytic(split, t, y, pfpt);
		break;
	case D_AUTODIFF:
		DtAutodiff(split, t, y, pfpt);
		break;
	case D_FORWARD:
		DtForward(split, t, y, pfpt);
		break;
	case D_CENTRED:
		DtCentred(split, t, y, pfpt);
		break;
	}
}

const BaseMat<FP>* BaseIVP::SplitMat(const FP t, const Vec<FP>& y, unsigned short split) {
	if( _jacSplitting ) {
		if( split == 1 )
			return Jac(t,y);
		throw Exception() << "Jacobian splitting does not define split matrix for split " << split << ".\n";
	}

	if( _splitMats[split] )
		delete _splitMats[split];
	_splitMats[split] = new Mat<FP>(y.Size(),y.Size());
	PhysicalSplitMat(split, t, y, (Mat<FP>&)*_splitMats[split]);
	return _splitMats[split];
}

const BaseMat<FP>* BaseIVP::Jac(const FP t, const Vec<FP>& y, unsigned short split) {
	if( _jacSplitting ) {
		if( split > 1 )
			throw Exception() << "Nonlinear term of Jacobian splitting has zero jacobian.";
		else
			split = 0;
	}

	if( !_jacFrozen ) {
		if( _splitJacs[split] )
			delete _splitJacs[split];
		_splitJacs[split] = new Mat<FP>(y.Size(), y.Size());

		switch( _jacType ) {
		case D_ANALYTIC:
			JacAnalytic(split, t, y, (Mat<FP>&)*_splitJacs[split]);
			break;
		case D_AUTODIFF:
			JacAutodiff(split, t, y, (Mat<FP>&)*_splitJacs[split]);
			break;
		case D_FORWARD:
			JacForward(split, t, y, (Mat<FP>&)*_splitJacs[split]);
			break;
		case D_CENTRED:
			JacCentred(split, t, y, (Mat<FP>&)*_splitJacs[split]);
			break;
		}
	}

	if( _jacScaling != 1. )
		*(Mat<FP>*)_splitJacs[split] *= _jacScaling;

	return _splitJacs[split];
}

const BaseMat<FP>* BaseIVP::SplitMatSparse(const FP t, const Vec<FP>& y, unsigned short split) {
	if( _jacSplitting ) {
		if( split == 1 )
			return JacSparse(t,y);
		throw Exception() << "Jacobian splitting does not define split matrix for split " << split << ".\n";
	}

	if( _splitMats[split] )
		delete _splitMats[split];
	_splitMats[split] = new CSRMat<FP>;
	PhysicalSplitMatSparse(split, t, y, (CSRMat<FP>&)*_splitMats[split]);
	return _splitMats[split];
}

const BaseMat<FP>* BaseIVP::JacSparse(const FP t, const Vec<FP>& y, unsigned short split) {
	if( _jacSplitting ) {
		if( split > 1 )
			throw Exception() << "Nonlinear term of Jacobian splitting has zero jacobian.";
		else
			split = 0;
	}

	if( !_jacFrozen ) {
		if( _splitJacs[split] )
			delete _splitJacs[split];
		_splitJacs[split] = new CSRMat<FP>;

		switch( _jacType ) {
		case D_ANALYTIC:
			JacAnalyticSparse(split, t, y, (CSRMat<FP>&)*_splitJacs[split]);
			break;
		case D_AUTODIFF:
			JacAutodiffSparse(split, t, y, (CSRMat<FP>&)*_splitJacs[split]);
			break;
		case D_FORWARD:
			JacForwardSparse(split, t, y, (CSRMat<FP>&)*_splitJacs[split]);
			break;
		case D_CENTRED:
			JacCentredSparse(split, t, y, (CSRMat<FP>&)*_splitJacs[split]);
			break;
		}
	}

	if( _jacScaling != 1. )
		*(CSRMat<FP>*)_splitJacs[split] *= _jacScaling;

/*	JacAutodiffSparse(split, t, y, (CSRMat<FP>&)*_splitJacs[split]);
	CSRMat<FP> ad((CSRMat<FP>&)*_splitJacs[split]);
	ad.ToDense().PrintMatlab(std::cout << "Autodiff = ", 8);

	JacAnalyticSparse(split, t, y, (CSRMat<FP>&)*_splitJacs[split]);
	CSRMat<FP> an((CSRMat<FP>&)*_splitJacs[split]);
	an.ToDense().PrintMatlab(std::cout << "Analytic = ", 8);

	std::cout << "error = " <<(ad-an).ToDense().MaxNorm() << "\n";
	throw Exception() << "Done.";*/

	return _splitJacs[split];
}

void BaseIVP::GetStats(Hash<ParamValue>& params) const {
}

void BaseIVP::PrintStats() const {
}

long BaseIVP::Size() const {
	return _initialCondition.Size();
}

#ifdef USE_ADOL_C
void BaseIVP::ToADOLC(const Vec<FP>& vec, Vec<adouble>& adolc) {
	for( long i = 0; i < vec.Size(); i++ )
		adolc[i] <<= vec[i];
}

void BaseIVP::FromADOLC(Vec<adouble>& adolc, Vec<FP>& vec) {
	for( long i = 0; i < vec.Size(); i++ )
		adolc[i] >>= vec[i];
}
#endif

