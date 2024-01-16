#ifndef _SINKHORN_SCALING_H
#define _SINKHORN_SCALING_H
#include "method.h"
/*
 * SinkhornScaling
 * Refernce
 * "Diagonal equivalence to matrices with prescribed row and column sums"
 * by R.Sinkhorn
 */
template<typename Solver>
class SinkhornScaling : public Methods
{
	CSRMatrix B;
	std::vector<double> L, R;
	mutable std::vector<double> Lb, y;
	Solver S;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","SinkhornScaling");
		ret.Set("maxiters",5);
		ret.Set("level","*");
		ret.Set("verbosity",1);
		ret.Set("write_matrix",0);
		ret.SubParameters("Solver") = Solver::DefaultParameters();
		return ret;
	}
	SinkhornScaling() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		B = A;
		bool print = GetParameters().template Get<int>("verbosity")? true : false;
		bool write_matrix = GetParameters().template Get<int>("write_matrix") ? true : false;
		int maxiters = GetParameters().template Get<int>("maxiters");
		int level = GetParameters().template Get<int>("level");
		idx_t size = A.Size();
		L.resize(size,1);
		R.resize(size,1);
		std::vector<double> u(size,0), v(size,0);
		double norm;
		for(int l = 0; l < maxiters; ++l)
		{
			if( print )
			{
				norm = 0;
				for(idx_t i = 0; i < size; ++i)
				{
					for(idx_t j = 0; j < B.RowSize(i); ++j)
					{
						double q = B.Val(i,j);
						q*= q;
						norm += q;
					}
				}
				std::cout << "mean norm: " << sqrt(norm/size) << std::endl;
			}
			for(idx_t i = 0; i < size; ++i)
			{
				u[i] = 0;
				for(idx_t j = 0; j < B.RowSize(i); ++j)
				{
					double q = B.Val(i,j);
					q*= q;
					u[i] += q;
				}
			}
			for(idx_t i = 0; i < size; ++i)
			{
				u[i] = 1.0/sqrt(u[i]);
				for(idx_t j = 0; j < B.RowSize(i); ++j)
					B.Val(i,j) *= u[i];
			}
			for(idx_t i = 0; i < size; ++i) v[i] = 0.0;
			for(idx_t i = 0; i < size; ++i)
			{
				for(idx_t j = 0; j < B.RowSize(i); ++j)
				{
					idx_t k = B.Col(i,j);
					double q = B.Val(i,j);
					q*= q;
					v[k] += q;
				}
			}
			for(idx_t i = 0; i < size; ++i) v[i] = 1.0/sqrt(v[i]);
			for(idx_t i = 0; i < size; ++i)
			{
				for(idx_t j = 0; j < B.RowSize(i); ++j)
					B.Val(i,j) *= v[B.Col(i,j)];
			}
			for(idx_t i = 0; i < size; ++i)
			{
				L[i] *= u[i];
				R[i] *= v[i];
			}
		}
		if( print )
		{
			norm = 0;
			for(idx_t i = 0; i < size; ++i)
			{
				for(idx_t j = 0; j < B.RowSize(i); ++j)
				{
					double q = B.Val(i,j);
					q*= q;
					norm += q;
				}
			}
			std::cout << "final mean norm: " << sqrt(norm/size) << std::endl;
		}
		if( write_matrix )
		{
			if( level > 0 )
			{
				if( print ) std::cout << "save B"<<level<<".mtx" << std::endl;
				B.Save("B"+to_string(level)+".mtx");
				if( print ) std::cout << "save L"<<level<<".mtx" << std::endl;
				SaveVector("L"+to_string(level)+".mtx",L);
				if( print ) std::cout << "save R"<<level<<".mtx" << std::endl;
				SaveVector("R"+to_string(level)+".mtx",R);
			}
			else 
			{
				if( print ) std::cout << "save B.mtx" << std::endl;
				B.Save("B.mtx");
				if( print ) std::cout << "save L.mtx" << std::endl;
				SaveVector("L.mtx",L);
				if( print ) std::cout << "save R.mtx" << std::endl;
				SaveVector("R.mtx",R);
			}
		}
		S.SetParameters(GetParameters().SubParameters("Solver"));
		bool success = S.Setup(B);
		Lb.resize(A.Size());
		y.resize(A.Size());
		return success;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		idx_t size = B.Size();
		// A x = b
		// L A x = L b
		// (L A R) y = L b, x = R y
		// B y = L b
		// y = B^{-1} L b
		// x = R y
		assert(b.size() == size);
		bool ret = false;
		x.resize(size,0.0);
		for(idx_t k = 0; k < size; ++k)
		{
			Lb[k] = L[k]*b[k];
			y[k] = x[k]/R[k];
		}
		ret = S.Solve(Lb,y);
		for(idx_t k = 0; k < size; ++k)
			x[k] = y[k]*R[k];
		return ret;
	}
	size_t Bytes() const {return S.Bytes() + B.Bytes() + get_bytes(Lb) + get_bytes(y);}
};

#endif //_SYMMETRICSCALING_H
