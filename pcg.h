#ifndef _PCG_H
#define _PCG_H
#include "method.h"
#include <cmath>

template<typename Preconditioner>
class PCG : public Methods
{
	Preconditioner P;
	mutable std::vector<double> r, p, z, w;
	const CSRMatrix * ptr_A;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","PCG");
		ret.Set("tol",1.0e-10);
		ret.Set("rtol",1.0e-7);
		ret.Set("dtol",1.0e+10);
		ret.Set("maxiters",5000);
		ret.Set("verbosity",1);
		ret.Set("true_residual", 1);
		ret.SubParameters("Preconditioner") = Preconditioner::DefaultParameters();
		return ret;
	}
	PCG() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		ptr_A = &A;
		P.SetParameters(GetParameters().SubParameters("Preconditioner"));
		r.resize(A.Size());
		p.resize(A.Size());
		z.resize(A.Size());
		w.resize(A.Size());
		return P.Setup(A);
	}
	bool ApplyPreconditioner(const std::vector<double> & r, std::vector<double> & z) const
	{
		Zero(z);
		if( !P.Solve(r,z) ) 
		{
			std::cout << "Preconditioner failed" << std::endl;
			return false;
		}
		return true;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		bool ptrue = GetParameters().Get<int>("true_residual") ? true : false;
		int maxiters = GetParameters().Get<int>("maxiters");
		double tol  = GetParameters().Get<double>("tol");
		double rtol = GetParameters().Get<double>("rtol");
		double dtol = GetParameters().Get<double>("dtol");
		int iters = 1;
		double resid, resid0, beta, alpha, kappa,ftol;
		const CSRMatrix & A = *ptr_A;
		idx_t size = A.Size();
        x.resize(size,0.0);
        std::copy(b.begin(),b.end(),r.begin());
        A.Multiply(-1.0,x,1.0,r);
        ApplyPreconditioner(r,z);
        std::copy(z.begin(),z.end(),p.begin()); //~ p = z;
		kappa = Dot(r,z);
		resid0 = resid = sqrt(fabs(kappa));
		ftol = std::max(tol,rtol*resid0);
		if (print)
		{
			std::cout << "PCG " << std::setw(4) << 0 << " " << std::setw(12) << resid << " | " << ftol;
			if( ptrue )
				std::cout << " true " << std::setw(12) << Resid(A, b, x);
			std::cout << std::endl;
			//~ std::cout << "\r";
			//~ std::cout.flush();
		}
        while( resid > ftol && resid < dtol && iters < maxiters+1 )
		{
			A.Multiply(1.0,p,0.0,w);//~ w = A*p;
			alpha = kappa/Dot(w,p);
			for(idx_t k = 0; k < size; ++k)
			{
				x[k] += alpha*p[k];
				r[k] -= alpha*w[k];
			}
			ApplyPreconditioner(r,z);
			beta = 1.0/kappa;
			kappa = Dot(r,z);
			beta *= kappa;
			resid = std::sqrt(fabs(kappa));
			for(idx_t k = 0; k < size; ++k)
				p[k] = z[k] + beta*p[k];
			if (print)
			{
				std::cout << "PCG "  << std::setw(4) << iters << " " << std::setw(14) << resid << " | " << ftol;
				if( ptrue )
					std::cout << " true " << std::setw(14) << Resid(A, b, x);
				std::cout << std::endl;
				//~ std::cout << "\r";
				//~ std::cout.flush();
			}
			iters++;
		}
		if (print)
		{
			std::cout << "PCG " << std::setw(4) << iters;
			if( ptrue )
				std::cout << " true " << std::setw(14) << Resid(A, b, x);
			std::cout << std::endl;
		}
		return resid <= ftol;
	}
	size_t Bytes() const {return sizeof(const CSRMatrix *) + P.Bytes() + get_bytes(r) + get_bytes(p) + get_bytes(z) + get_bytes(w);}
};

#endif //_PCG_H
