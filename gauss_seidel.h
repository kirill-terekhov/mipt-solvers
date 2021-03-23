#ifndef _GAUSS_SEIDEL_H
#define _GAUSS_SEIDEL_H
#include "method.h"
/*
 * GaussSeidel
 * this is symmetric gauss-seidel method if number of iterations is even,
 * or symmetric successive over relaxation if omega is not unit.
 */

class GaussSeidel : public Methods
{
	std::vector<idx_t> d;
	const CSRMatrix * ptr_A;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","GaussSeidel");
		ret.Set("tol",0.0);
		ret.Set("omega",1.0);
		ret.Set("maxiters",2);
		ret.Set("verbosity",0);
		return ret;
	}
	GaussSeidel() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		ptr_A = &A;
		//search diagonal
		if( !A.DiagonalPosition(d) )
		{
			std::cout << "No diagonal elements in matrix" << std::endl;
			return false;
		}
		return true;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool print   = GetParameters().Get<int>("verbosity") & 2 ? true : false;
		int maxiters = GetParameters().Get<int>("maxiters");
		double tol   = GetParameters().Get<double>("tol");
		double omega = GetParameters().Get<double>("omega");
		const CSRMatrix & A = *ptr_A;
		double resid = 0;
		int iters = 0;
		x.resize(A.Size(),0.0);
		do
		{
			for(idx_t it = 0; it < A.Size(); ++it)
			{
				//idx_t i = it; //forward substitution only
				//idx_t i = A.Size()-1-it; //backward substitution only
				idx_t i = iters%2 ? A.Size()-1-it : it; //alternate backward and forward subsitution
				double s = 0.0;
				for(idx_t j = 0; j < d[i]; ++j)
					s += A.Val(i,j)*x[A.Col(i,j)];
				for(idx_t j = d[i]+1; j < A.RowSize(i); ++j)
					s += A.Val(i,j)*x[A.Col(i,j)];
				x[i] = omega*(b[i] - s)/A.Val(i,d[i]) + (1.0-omega)*x[i];
			}
			if( tol )
			{
				resid = Resid(A,b,x);
				if( resid < tol ) break;
			}
			iters++;
		} while( iters < maxiters );
		if( print )
			std::cout << "Gauss-Seidel " << iters << " true " << (tol ? resid : Resid(A,b,x)) << std::endl;
		return !tol || resid < tol;
	}
	size_t Bytes() const {return sizeof(const CSRMatrix *) + get_bytes(d);}
};

#endif //_GAUSS_SEIDEL_H
