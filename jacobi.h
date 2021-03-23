#ifndef _JACOBI_H
#define _JACOBI_H
#include "method.h"

/* Jacobi 
 * this is jacobi method with relaxation parameter.
 */

class Jacobi : public Methods
{
	std::vector<idx_t> d;
	mutable std::vector<double> x0;
	const CSRMatrix * ptr_A;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","Jacobi");
		ret.Set("tol",0.0);
		ret.Set("maxiters",2);
		ret.Set("verbosity",0);
		ret.Set("omega",1.0);
		return ret;
	}
	Jacobi() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		idx_t err = 0;
		ptr_A = &A;
		//search diagonal
		if( !A.DiagonalPosition(d) )
		{
			std::cout << "No diagonal elements in matrix" << std::endl;
			return false;
		}
		x0.resize(A.Size());
		return err == 0;
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
			std::copy(x.begin(),x.end(),x0.begin());
			for(idx_t i = 0; i < A.Size(); ++i)
			{
				double s = 0.0;
				for(idx_t j = 0; j < d[i]; ++j)
					s += A.Val(i,j)*x0[A.Col(i,j)];
				for(idx_t j = d[i]+1; j < A.RowSize(i); ++j)
					s += A.Val(i,j)*x0[A.Col(i,j)];
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
			std::cout << "Jacobi " << iters << " true " << (tol? resid : Resid(A,b,x)) << std::endl;
		return !tol || resid < tol;
	}
	size_t Bytes() const {return sizeof(const CSRMatrix *) + get_bytes(x0) + get_bytes(d);}
};

#endif //_JACOBI_H
