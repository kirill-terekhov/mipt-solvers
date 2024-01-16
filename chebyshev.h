#ifndef _CHEBYSHEV_H
#define _CHEBYSHEV_H
#include "method.h"
#include <cmath>

/* Chebyshev 
 * this is chebyshev polynomial smoother with
 * the eigenvalues estimated using Gershgorin's disks.
 * 
 * [1] Improving the arithmetic intensity of multigrid with the help of polynomial smoothers
 * by P. Ghysels, P. Klosiewicz, W. Vanroose
 * [2] Parallel multigrid smoothing: polynomial versus Gauss–Seidel
 * by M. Adams, M. Brezina, J. Hu, R. Tuminaro
 */

class Chebyshev : public Methods
{
	mutable std::vector<double> p;
	const CSRMatrix * ptr_A;
	double d, c;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","Chebyshev");
		ret.Set("tol",0.0);
		ret.Set("maxiters",2);
		ret.Set("verbosity",1);
		return ret;
	}
	Chebyshev() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		bool print   = GetParameters().Get<int>("verbosity") & 1 ? true : false;
		idx_t err = 0;
		ptr_A = &A;
		double lmin = 1.0e+20, lmax = -1.0e+20;
		//estimate large / small eigenvalue
		{ // Girshgorin circle theorem
			std::vector<double> D(A.Size(),0.0);
			//get out the diagonal
			A.Diagonal(D);
			//compute the radius and estimate eigenvalues
			for(idx_t k = 0; k < A.Size(); ++k)
			{
				double R = 0;
				for(idx_t j = 0; j < A.RowSize(k); ++j) if( A.Col(k,j) != k )
					R += std::abs(D[A.Col(k,j)]*A.Val(k,j)/D[k]);
					//~ R += std::abs(A.Val(k,j));
				lmin = std::min(lmin,D[k] - R);
				lmax = std::max(lmax,D[k] + R);
			}
			if( std::abs(lmax) < std::abs(lmin) )
				std::swap(lmax,lmin);
			if( print )
			{
				std::cout << "Estimated bounds of eigenvalues: " << lmin << ":" << lmax << "." << std::endl;
				std::cout << "Consumed by D: " << get_bytes(D)/1024.0/1024.0 << "Mb." << std::endl;
			}
		}
		lmin = lmax/30.0;
		d = 0.5*(lmax+lmin);
		c = 0.5*(lmax-lmin);
		if( print )
			std::cout << "Eigenvalues: " << lmin << ":" << lmax << " center: " << d << " foci:" << c << std::endl;
		p.resize(A.Size());
		return err == 0;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool print   = GetParameters().Get<int>("verbosity") & 2 ? true : false;
		int maxiters = GetParameters().Get<int>("maxiters");
		double tol   = GetParameters().Get<double>("tol");
		const CSRMatrix & A = *ptr_A;
		idx_t size = (idx_t)A.Size();
		double resid = 0, alpha, beta;
		int iters = 0;
		x.resize(A.Size(),0.0);
		std::fill(p.begin(),p.end(),0.0);
		do
		{
			if( iters == 0 )
				alpha = 1.0/d;
			else if( iters == 1 )
				alpha = 2*d / (2*d*d - c*c);
			else
				alpha = 1.0/(d - 0.25 * alpha * c * c);
			beta = alpha*d - 1.0;
			for(idx_t k = 0; k < size; ++k)
			{
				double rk = b[k];
				for(idx_t j = 0; j < A.RowSize(k); ++j)
					rk -= A.Val(k,j)*x[A.Col(k,j)];
				p[k] = alpha*rk + beta*p[k];
			}
			for(idx_t k = 0; k < size; ++k)
				x[k] += p[k];
			if( tol )
			{
				resid = Resid(A,b,x);
				if( print )
				{
					std::cout << "Chebyshev " << iters << " true " << resid;
					std::cout << "/r";
					//~ std::cout << std::endl;
				}
				if( resid < tol ) break;
			}
			iters++;
		} while( iters < maxiters );
		if( print )
			std::cout << "Chebyshev " << iters << " true " << (tol? resid : Resid(A,b,x)) << std::endl;
		return !tol || resid < tol;
	}
	size_t Bytes() const {return sizeof(const CSRMatrix *) + get_bytes(p) + 2*sizeof(double);}
};

#endif //_CHEBYSHEV_H
