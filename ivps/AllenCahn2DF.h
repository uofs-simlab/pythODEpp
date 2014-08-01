#ifndef ALLEN_CAHN_2DF_H
#define ALLEN_CAHN_2DF_H

class AllenCahn2DF : public TwoSplittingIVP
{
protected:
	FP alpha;
	FP gamma;
	long N;
	FP delta;
	FP d;

	void JacAnalyticSparse(unsigned short split, const FP t, const Vec<FP>& y,
		CSRMat<FP>& jac)
	{
		if(split > 2)
		{
			throw Exception() << "Analytic split " << split << " is not defined.";
		}

		FP* values;
		long* columns;
		long* rowIndex;
		long N2 = N * N;
		FP tmp;

		if(split == 2)
		{			
			values = new FP[N2];
			columns = new long[N2];
			rowIndex = new long[N2];

			for(long i = 0; i < N2; i++)
			{
				tmp = y[i];
				values[i] = gamma * (1 - 3 * tmp * tmp);
				columns[i] = i;
				rowIndex[i] = i;
			}

			jac = CSRMat<FP>(values, columns, rowIndex, N2, N2, N2);
		}
		else
		{
			long nnz = 5 * N2 - 4 * N;
			values = new FP[nnz];
			columns = new long[nnz];
			rowIndex = new long[N2];
			rowIndex[0] = 0;

			long ind = 0;
			for(long i = 0; i < N2; i++)
			{
				if(i >= N)
				{
					values[ind] = d;
					columns[ind] = i - N;
					rowIndex[i] = ind;
					ind = ind + 1;
				}
				if(i % N != 0)
				{
					values[ind] = d;
					columns[ind] = i - 1;
					if(i < N)
						rowIndex[i] = ind;
					ind = ind + 1;
				}
				
				values[ind] = -4 * d;
				if(split == 0)
				{
					tmp = y[i];
					values[ind] = values[ind] + gamma * (1 - 3 * tmp * tmp);
				}
				columns[ind] = i;
				ind = ind + 1;

				if(i % N != N - 1)
				{
					values[ind] = d;
					columns[ind] = i + 1;
					ind = ind + 1;
				}
				if(i < N2 - N)
				{
					values[ind] = d;
					columns[ind] = i + N;
					ind = ind + 1;
				}
			}

			jac = CSRMat<FP>(values, columns, rowIndex, N2, N2, nnz);
		}

		delete[] values;
		delete[] columns;
		delete[] rowIndex;
	}

	template <class T>
	void Split1Internal(const T t, const Vec<T>& y, Vec<T>& yp)
	{
		T l, r, u, d;
		long index;

		for(long i = 0; i < N; i++)
		{
			for(long j = 0; j < N; j++)
			{
				index = N * i + j;
				
				if(index <= N - 1)
					d = exactSolution((T) (delta * (j + 1)), (T) 0, t);
				else
					d = y[index - N];

				if(index % N == 0)
					l = exactSolution((T) 0, (T) (delta * (i + 1)), t);
				else
					l = y[index - 1];

				if(index % N == N - 1)
					r = exactSolution((T) 1, (T) (delta * (i + 1)), t);
				else
					r = y[index + 1];

				if(index >= N * N - N)
					u = exactSolution((T) (delta * (j + 1)), (T) 1, t);
				else
					u = y[index + N];

				yp[index] = (u + d + l + r - 4 * y[index]) * alpha / delta
					/ delta;
			}
		}
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp)
	{
		long index;
		T yind;

		for(long i = 0; i < N; i++)
		{
			for(long j = 0; j < N; j++)
			{
				index = N * i + j;
				yind = y[index];
				yp[index] = gamma * yind * (1 - yind * yind)
					+ f(delta * (j + 1), delta * (i + 1), t);
			}
		}
	}

	template <class T>
	T f(const T x, const T y, const T t)
	{
		T xt = M_PI * 2 * (x - t);
		T yt = M_PI * 3 * (y - t);
		T sxt = sin(xt);
		T cxt = cos(xt);
		T syt = sin(yt);
		T cyt = cos(yt);
		
		T result = M_PI * (3 * sxt * syt - 2 * cxt * cyt);
		result = result + gamma * sxt * sxt * cyt * cyt * (sxt * cyt + 6);
		result = result + sxt * cyt * (11 * gamma + 13 * alpha * M_PI * M_PI);
		result = result + 6 * gamma;

		return result;
	}

	template <class T>
	T exactSolution(const T x, const T y, const T t)
	{
		return sin(M_PI * 2 * (x - t)) * cos(M_PI * 3 * (y - t)) + 2;
	}

public:
	AllenCahn2DF(Hash<ParamValue>& params) : TwoSplittingIVP(params)
	{
		alpha = GetDefaultFP(params, "alpha", 0.01);
		gamma = GetDefaultFP(params, "gamma", 3);
		N = GetDefaultLong(params, "N", 25);
		delta = 1. / (N + 1);
		d = alpha * (N + 1) * (N + 1);

		_initialCondition.Resize(N * N);
		for(long i = 0; i < N; i++)
		{
			for(long j = 0; j < N; j++)
			{
				_initialCondition[N * i + j] = exactSolution(delta * (j + 1),
					delta * (i + 1), 0.);
			}
		}
	}

	LINK_TWOSPLIT
	IVP_NAME("Allen Cahn 2DF")
};

#endif

