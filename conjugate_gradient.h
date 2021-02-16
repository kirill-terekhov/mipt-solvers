#ifndef _CONJUGATE_GRADIENT_H
#define _CONJUGATE_GRADIENT_H
#include "method.h"
#include <cmath>

class ConjugateGradient : public Methods
{
	const CSRMatrix * ptr_A;
	mutable std::vector<double> r, p, w;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","ConjugateGradient");
		ret.Set("tol",1.0e-7);
		ret.Set("maxiters",5000);
		ret.Set("verbosity",1);
		return ret;
	}
	ConjugateGradient() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		ptr_A = &A;
		return true; // test A is symmetric?
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool print = GetParameters().Get<int>("verbosity") & 2 ? true : false;
		int iters = 0, maxiters = GetParameters().Get<int>("maxiters");
		double resid, beta, alpha, tol = GetParameters().Get<double>("tol");
		const CSRMatrix & A = *ptr_A;
        x.resize(A.Size(),0.0);
        std::copy(b.begin(),b.end(),r.begin());
        A.Multiply(-1.0,x,0.0,r); //~ r -= A*x;
        std::copy(r.begin(),r.end(),p.begin()); //~ p = r;
        resid = Dot(r,r);
        while( sqrt(fabs(resid)) > tol && iters < maxiters )
		{
			A.Multiply(1.0,p,0.0,w);//~ w = A*p;
			alpha = resid/Dot(w,p);
			for(idx_t k = 0; k < A.Size(); ++k)
			{
				x[k] += alpha*p[k];
				r[k] -= alpha*w[k];
			}
			beta = 1.0/resid;
			resid = Dot(r,r);
			beta *= resid;
			for(idx_t k = 0; k < A.Size(); ++k)
				p[k] = r[k] + beta*p[k];
			if( print )
				std::cout << "CG " << iters << " " << sqrt(fabs(resid)) << " true " << Resid(A,b,x) << std::endl;
			iters++;
		}
		if( print ) 
			std::cout << "CG " << iters << " true " << Resid(A,b,x) << std::endl;
		return resid < tol;
	}
	size_t Bytes() const {return sizeof(const CSRMatrix *) + get_bytes(r) + get_bytes(p) + get_bytes(w);}
};

#endif //_CONJUGATE_GRADIENT_H
