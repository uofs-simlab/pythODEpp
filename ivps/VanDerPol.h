#ifndef VAN_DER_POL_H
#define VAN_DER_POL_H

class VanDerPol : public TwoSplittingIVP
{
protected:
	FP epsilon;

	void JacAnalytic(unsigned short split, const FP t, const Vec<FP>& y, Mat<FP>& jac)
	{
		FP Y = y[0];
		FP Z = y[1];

		jac(0, 0) = 0;

		if(split == 0 || split == 1)
		{
			jac(1, 0) = (-2 * Y * Z - 1) / epsilon;
			jac(1, 1) = (1 - Y * Y) / epsilon;
			if(split == 1)
			{
				jac(0, 1) = 0;
			}
		}
		
		if(split == 0 || split == 2)
		{
			jac(0, 1) = 1;
			if(split == 2)
			{
				jac(1, 0) = 0;
				jac(1, 1) = 0;
			}
		}

		if(split > 2)
		{
			throw Exception() << "Analytic split " << split << " is not defined.";
		}
	}

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp)
	{
		T Y = y[0];
		T Z = y[1];
		yp[0] = 0;
		yp[1] = ((1 - Y * Y) * Z - Y) / epsilon;
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp)
	{
		yp[0] = y[1];
		yp[1] = 0;
	}

public:
	VanDerPol(Hash<ParamValue>& params) : TwoSplittingIVP(params)
	{
		epsilon = GetDefaultFP(params, "epsilon", 0.01);
		_initialCondition.Resize(2);
		_initialCondition[0] = GetDefaultFP(params, "y0", 2);
		_initialCondition[1] = GetDefaultFP(params, "z0", -0.6654321);
	}

	LINK_TWOSPLIT
	IVP_NAME("Van der Pol Equation")
};

#endif

