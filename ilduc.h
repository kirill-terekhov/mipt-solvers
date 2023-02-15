#ifndef _ILDUC_H
#define _ILDUC_H


#include "csrmatrix.h"
#include "csrtriangular.h"
#include "csctraversal.h"
#include "row_accumulator.h"
#include "invnorm.h"
#include "method.h"

/*
 * ILDUC - Crout incomplete LDU factorization
 * Implemented following 
 * [1] "Crout Versions of ILU for General Sparse Matrices"
 * by N.Li, Y.Saad, E.Show
 * [2] "A robust and efficient ILU that incorporates the growth of the inverse triangular factors"
 * by M.Bollhoefer
 */
 
 
class ILDUC : public Methods
{
	CSRTriangular L, U;
	std::vector<double> D;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","ILDUC");
		ret.Set("drop_tolerance",0.01);
		ret.Set("diagonal_tolerance",1.0e-18);
		ret.Set("diagonal_perturbation",0.0);
		ret.Set("write_matrix",0);
		ret.Set("verbosity",1);
		ret.Set("inverse_estimation",1);
		ret.Set("premature_dropping",0);
		ret.Set("check",0);
		return ret;
	}
	ILDUC() : L(CSRTriangular::LowerCSC), U(CSRTriangular::UpperCSR) {GetParameters() = DefaultParameters();}
	~ILDUC() {}
	bool Setup(const CSRMatrix & Ain)
	{
		bool print = GetParameters().Get<int>("verbosity") & 1 ? true : false;
		bool check = GetParameters().Get<int>("check") ? true : false;
		bool invest = GetParameters().Get<int>("inverse_estimation") ? true : false;
		bool predrop = GetParameters().Get<int>("premature_dropping") ? true : false;
		bool write_matrix = GetParameters().Get<int>("write_matrix") ? true : false;
		double tau = GetParameters().Get<double>("drop_tolerance");
		double pert = GetParameters().Get<double>("diagonal_perturbation");
		double dtol = GetParameters().Get<double>("diagonal_tolerance");
		idx_t report_pace = std::max<idx_t>(1,Ain.Size() / 25);
		CSRMatrix A = Ain;
		L.Clear();
		U.Clear();
		A.SortRows();
		CSCTraversal At(A,A.Size()), Lt(L,A.Size()), Ut(U,A.Size());
		RowAccumulator<double> u(A.Size()), l(A.Size());
		std::vector<double> unorms(A.Size(),0.0), lnorms(A.Size(),0.0);
		Invnorm iLest(invest | print ? A.Size() : 0), iUest(invest | print ? A.Size() : 0);
		double iLnorm = 1, iUnorm = 1, Dmax = 0, Dmin = 1.0e+20;
		D.resize(A.Size(),0.0);
		if( print ) 
			std::cout << "tau: " << tau << " inverse estimation: " << (invest ? "yes":"no") << std::endl;
		for (idx_t k = 0; k < A.Size(); k++)
		{
			if( print )
			{
				if( k % report_pace == 0 )
				{
					std::cout << "precond: " << std::setw(12) << (k+1)*100.0 / A.Size() << " L " << iLnorm << " D " << Dmax/Dmin << " U " << iUnorm << "\r";
					std::cout.flush();
				}
			}
			{
				//Compute U-part
				{
					//Uncompress row
					idx_t curr = k;
					u[k] = 0.0;
					for(idx_t j = 0; j < A.RowSize(k); ++j) if( A.Col(k,j) >= k )
						curr = u.InsertOrdered(curr,A.Col(k,j),A.Val(k,j));
					if( pert )
						u.Get(k) = u.Get(k) * (1.0 + pert) + (u.Get(k) < 0.0 ? -1.0 : 1.0) * pert;
					//U part elimination with L
					for(idx_t j = Lt.Begin(); j != Lt.End(); j = Lt.Next(j))
					{
						double v = -L.Val(j,Lt.Position(j)) * D[j];
						if( !predrop || std::fabs(v)*unorms[j] > tau )
						{
							curr = k;
							for(idx_t r = 0; r < U.RowSize(j); ++r) if( U.Col(j,r) >= k )
								curr = u.InsertOrdered(curr,U.Col(j,r),v*U.Val(j,r));
						}
					}
				}
				//Compute L-part
				{
					//Uncompress column
					idx_t curr = k;
					l[k] = 0.0;
					for(idx_t j = At.Begin(); j != At.End(); j = At.Next(j)) if( j >= k )
						curr = l.InsertOrdered(curr,j,A.Val(j,At.Position(j)));
					At.NextColumn();
					if( pert )
						l.Get(k) = l.Get(k) * (1.0 + pert) + (l.Get(k) < 0.0 ? -1.0 : 1.0) * pert;
					//L part elimination with U
					for(idx_t j = Ut.Begin(); j != Ut.End(); j = Ut.Next(j))
					{
						double v = -U.Val(j,Ut.Position(j)) * D[j];
						if( !predrop || std::fabs(v)*lnorms[j] > tau )
						{
							curr = k;
							for(idx_t r = 0; r < L.RowSize(j); ++r) if( L.Col(j,r) >= k )
								curr = l.InsertOrdered(curr,L.Col(j,r),v*L.Val(j,r));
						}
					}
				}
			}
			//retrive diagonal
			D[k] = (u.Get(k) + l.Get(k))*0.5;
			D[k] = (D[k] < 0.0 ? -1.0 : 1.0)*std::max(std::fabs(D[k]),dtol);
			Dmax = std::max(Dmax,std::fabs(D[k]));
			Dmin = std::min(Dmin,std::fabs(D[k]));
			{
				//Assemble U-part
				{
					//Scale u by diagonal
					double unorm = 0;
					for(idx_t j = u.Begin(); j != u.End(); j = u.Next(j))
					{
						u.Get(j) /= D[k];
						unorm += std::pow(u.Get(j),2);
					}
					unorm = std::sqrt(unorm);
					unorms[k] = unorm;
					//inverse estimation
					if( invest | print )
					{
						iUnorm = iUest.Estimate(u,k);
						iUest.Update(u,k);
						if( invest ) unorm /= iUnorm;
					}
					//assemble row of U
					U.PushBack(k,u.Get(k));
					for(idx_t j = u.Next(u.Begin()); j != u.End(); j = u.Next(j))
					{
						double v = u.Get(j);
						if( std::fabs(v) > tau*unorm ) U.PushBack(j,v);
					}
					U.FinalizeRow();
					Ut.NewRow(k);
					Ut.NextColumn();
					u.Clear();
				}
				//Assemble L-part
				{
					//Scale u by diagonal
					double lnorm = 0;
					for(idx_t j = l.Begin(); j != l.End(); j = l.Next(j))
					{
						l.Get(j) /= D[k];
						lnorm += std::pow(l.Get(j),2);
					}
					lnorm = std::sqrt(lnorm);
					lnorms[k] = lnorm;
					//inverse estimation
					if( invest | print )
					{
						iLnorm = iLest.Estimate(l,k);
						iLest.Update(l,k);
						if( invest ) lnorm /= iLnorm;
					}
					//assemble row of U
					L.PushBack(k,l.Get(k));
					for(idx_t j = l.Next(l.Begin()); j != l.End(); j = l.Next(j))
					{
						double v = l.Get(j);
						if( std::fabs(v) > tau*lnorm ) L.PushBack(j,v);
					}
					L.FinalizeRow();
					Lt.NewRow(k);
					Lt.NextColumn();
					l.Clear();
				}
			}
		}
		if( print )
		{
			std::cout << "      nonzeros in A: " << A.Nonzeros() << " consumed: " << A.Bytes()/1024.0/1024.0 << "Mb " << std::endl;
			std::cout << "      nonzeros in L: " << L.Nonzeros() << " consumed: " << L.Bytes()/1024.0/1024.0 << "Mb " << std::endl;
			std::cout << "      nonzeros in U: " << U.Nonzeros() << " consumed: " << U.Bytes()/1024.0/1024.0 << "Mb " << std::endl;
			std::cout << "      consumed by D: " << D.capacity()*sizeof(double)/1024.0/1024.0 << "Mb " << std::endl;
			size_t bytes = A.Bytes() + L.Bytes() + U.Bytes();
			bytes += (unorms.capacity()+lnorms.capacity()+D.capacity())*sizeof(double);
			bytes += Lt.Bytes() + Ut.Bytes();
			bytes += iUest.Bytes() + iLest.Bytes();
			std::cout << "      consumed in factorization " << bytes/1024.0/1024.0 << "Mb " << std::endl;
			std::cout << "      fill-in LU: " << (L.Nonzeros() + U.Nonzeros())/(1.0*A.Nonzeros()) << std::endl;
			std::cout << "      estimated inverse norms L " << iLnorm << " D " << Dmax/Dmin << " U " << iUnorm << std::endl;
		}
		if( write_matrix )
		{
			if( print ) std::cout << "save L.mtx" << std::endl;
			L.Save("L.mtx");
			if( print ) std::cout << "save U.mtx" << std::endl;
			U.Save("U.mtx");
		}
		if( check )
		{
			if( !L.Check() )
				std::cout << "L is not lower-triangular" << std::endl;
			if( !U.Check() )
				std::cout << "U is not upper-triangular" << std::endl;
			CSRMatrix B = L.Transpose() * CSRMatrix(D) * U;
			double difmax = -1.0e+20, difmin = 1.0e+20;
			for (idx_t k = 0; k < A.Size(); ++k)
			{
				double Asum = 0, Bsum = 0, sum = 0;
				for (idx_t j = 0; j < A.RowSize(k); ++j)
					Asum += A.Val(k, j);
				for (idx_t j = 0; j < B.RowSize(k); ++j)
					Bsum += B.Val(k, j);
				sum = Asum - Bsum;
				//if (std::fabs(sum) > 1.0e-1)
				//	std::cout << k << " row sum difference " << sum << " in A " << Asum << " in B " << Bsum << std::endl;
				difmax = std::max(sum, difmax);
				difmin = std::min(sum, difmin);
			}
			std::cout << "rowsum difference LU-A interval " << difmin << ":" << difmax << std::endl;
			CSRMatrix C = B - A;
			double err = 0, errmax = -1.0e+20, errmin = 1.0e+20;
			for (idx_t k = 0; k < C.Size(); ++k)
			{
				for (idx_t j = 0; j < C.RowSize(k); ++j)
				{
					errmax = std::max(errmax, C.Val(k, j));
					errmin = std::min(errmin, C.Val(k, j));
					err += pow(C.Val(k, j), 2);
				}
			}
			err = std::sqrt(err);
			
			std::cout << "matrix difference LU-A L2-norm " << err << " per-entry interval " << errmin << ":" << errmax << std::endl;
		}
		return true;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool print = GetParameters().Get<int>("verbosity") & 2 ? true : false;
		x.resize(b.size());
		std::copy(b.begin(),b.end(),x.begin());
		if( print ) std::cout << "ILDUC Solve with L" << std::endl;
		L.Solve(x);
		if( print ) std::cout << "ILDUC Solve with D" << std::endl;
		for(idx_t k = 0; k < D.size(); ++k) 
			x[k] /= D[k];
		if( print ) std::cout << "ILDUC Solve with U" << std::endl;
		U.Solve(x);
		return true;
	}
	size_t Bytes() const {return U.Bytes() + L.Bytes() + get_bytes(D) + sizeof(const CSRMatrix *);}
};

#endif //_ILDUC_H
