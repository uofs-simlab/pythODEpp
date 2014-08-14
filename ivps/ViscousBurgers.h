#ifndef VISCOUS_BURGERS_H
#define VISCOUS_BURGERS_H

class ViscousBurgers : public TwoSplittingIVP
{
protected:
	FP nu;
	long N;
	FP delta;

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
		long nnz = 5 * N2 - 4 * N;
		values = new FP[nnz];
		columns = new long[nnz];
		rowIndex = new long[N2];

		FP ndd = 0;
		if(split == 0 || split == 1)
			ndd = nu / delta / delta;

		FP l, r, u, d;
		rowIndex[0] = 0;
		long ind = 0;
		for(long i = 0; i < N2; i++)
		{
			if(i >= N)
			{
				values[ind] = ndd;
				if(split == 0 || split == 2)
					values[ind] = values[ind] + y[i] / 2 / delta;

				columns[ind] = i - N;
				rowIndex[i] = ind;
				ind = ind + 1;
			}
			if(i % N != 0)
			{
				values[ind] = ndd;
				if(split == 0 || split == 2)
					values[ind] = values[ind] + y[i] / 2 / delta;

				columns[ind] = i - 1;
				if(i < N)
					rowIndex[i] = ind;
				ind = ind + 1;
			}

			values[ind] = -4 * ndd;
			if(split == 0 || split == 2)
			{
				if(i <= N - 1)
					d = exactSolution(delta * (i + 1), 0., t);
				else
					d = y[i - N];

				if(i % N == 0)
					l = exactSolution(0., delta * (i / N + 1), t);
				else
					l = y[i - 1];

				if(i % N == N - 1)
					r = exactSolution(0.5, delta * (i / N + 1), t);
				else
					r = y[i + 1];

				if(i >= N * N - N)
					u = exactSolution(delta * (i % N + 1), 0.5, t);
				else
					u = y[i + N];

				values[ind] = values[ind] - (r - l + u - d) / 2 / delta;
			}
			columns[ind] = i;
			ind = ind + 1;

			if(i % N != N - 1)
			{
				values[ind] = ndd;
				if(split == 0 || split == 2)
					values[ind] = values[ind] - y[i] / 2 / delta;

				columns[ind] = i + 1;
				ind = ind + 1;
			}
			if(i < N2 - N)
			{
				values[ind] = ndd;
				if(split == 0 || split == 2)
					values[ind] = values[ind] - y[i] / 2 / delta;

				columns[ind] = i + N;
				ind = ind + 1;
			}
		}

		jac = CSRMat<FP>(values, columns, rowIndex, N2, N2, nnz);

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
					r = exactSolution((T) 0.5, (T) (delta * (i + 1)), t);
				else
					r = y[index + 1];

				if(index >= N * N - N)
					u = exactSolution((T) (delta * (j + 1)), (T) 0.5, t);
				else
					u = y[index + N];

				yp[index] = (u + d + l + r - 4 * y[index]) * nu / delta
					/ delta;
			}
		}
	}

	template <class T>
	void Split2Internal(const T t, const Vec<T>& y, Vec<T>& yp)
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
					r = exactSolution((T) 0.5, (T) (delta * (i + 1)), t);
				else
					r = y[index + 1];

				if(index >= N * N - N)
					u = exactSolution((T) (delta * (j + 1)), (T) 0.5, t);
				else
					u = y[index + N];

				yp[index] = -y[index] * (r - l + u - d) / 2 / delta;
			}
		}
	}

	template <class T>
	T exactSolution(const T x, const T y, const T t)
	{
		return 1.0 / (1 + exp((x + y - t) / 2 / nu));
	}

public:
	ViscousBurgers(Hash<ParamValue>& params) : TwoSplittingIVP(params)
	{
		nu = GetDefaultFP(params, "nu", 0.1);
		N = GetDefaultLong(params, "N", 25);
		delta = 0.5 / (N + 1);

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
	IVP_NAME("Viscous Burgers'")
};

#endif

