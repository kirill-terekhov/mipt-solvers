#ifndef _BICGSTAB_H
#define _BICGSTAB_H
#include "method.h"
#include <cmath>

template<typename Preconditioner>
class BICGSTAB : public Methods
{
	Preconditioner P;
	mutable std::vector<double> r, r0, y, z, p, v, s, t;
	const CSRMatrix * ptr_A;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","BICGSTAB");
		ret.Set("tol",1.0e-8);
		ret.Set("rtol",1.0e-6);
		ret.Set("dtol",1.0e+50);
		ret.Set("maxiters",5000);
		ret.Set("verbosity",1);
		ret.Set("true_residual", 0);
		ret.SubParameters("Preconditioner") = Preconditioner::DefaultParameters();
		return ret;
	}
	BICGSTAB() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		ptr_A = &A;
		P.SetParameters(GetParameters().SubParameters("Preconditioner"));
		r.resize(A.Size());
		r0.resize(A.Size());
		y.resize(A.Size());
		z.resize(A.Size());
		p.resize(A.Size());
		v.resize(A.Size());
		s.resize(A.Size());
		t.resize(A.Size());
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
		bool print = GetParameters().template Get<int>("verbosity")? true : false;
		bool ptrue = GetParameters().template Get<int>("true_residual") ? true : false;
		int maxiters = GetParameters().template Get<int>("maxiters");
		double tol  = GetParameters().template Get<double>("tol");
		double rtol = GetParameters().template Get<double>("rtol");
		double dtol = GetParameters().template Get<double>("dtol");
		int iters = 1;
		idx_t size = ptr_A->Size();
		const CSRMatrix & A = *ptr_A;
		double alpha = 1, beta = 0, omega = 1, rho = 1, resid, resid0, ftol;
		x.resize(size,0.0);
		std::fill(p.begin(),p.end(),0.0);
		std::fill(v.begin(),v.end(),0.0);
		std::copy(b.begin(),b.end(),r.begin());
		A.Multiply(-1.0,x,1.0,r); //~ r -= A*x;
		resid0 = resid = Norm(r);
		ftol = std::max(tol,rtol*resid0);
		if (print)
		{
			std::cout << "BICGSTAB " << std::setw(4) << 0 << " " << std::setw(12) << resid << " | " << ftol;
			if( ptrue )
				std::cout << " true " << std::setw(12) << Resid(A, b, x);
			//std::cout << "\r";
			std::cout << std::endl;
			std::cout.flush();
		}
		if( resid <= ftol )
			return true;
		std::copy(r.begin(),r.end(),r0.begin());
		do 
		{
			beta = 1.0 / rho * alpha / omega;
			rho = Dot(r0,r);
			if (std::fabs(rho) < 1.0e-15)
			{
				for (idx_t i = 0; i < size; ++i)
					r0[i] += (2.0 * rand() / (1.0 * RAND_MAX) - 1.0) / (1.0 * r0.size());
				resid0 = Norm(r0);
				for (idx_t i = 0; i < size; ++i)
					r0[i] /= resid0;
				rho = Dot(r0, r);
			}
			beta*= rho;
			for (idx_t i = 0; i < size; ++i)
				p[i] = r[i] + beta*(p[i] - omega*v[i]);
			ApplyPreconditioner(p,y);
			A.Multiply(1.0,y,0.0,v);//~ v = A*y;
			alpha = rho / Dot(r0,v);
			for (idx_t i = 0; i < size; ++i)
				s[i] = r[i] - alpha*v[i];
			ApplyPreconditioner(s,z);
			A.Multiply(1.0,z,0.0,t);//~ t = A*z;
			omega = Dot(t,s);
			if( omega )
			{
				omega /= Dot(t,t);
				for (idx_t i = 0; i < size; ++i)
				{
					x[i] += alpha*y[i] + omega*z[i];
					r[i] = s[i] - omega*t[i];
				}
			}
			else
			{
				for (idx_t i = 0; i < size; ++i)
				{
					x[i] += alpha*y[i];
					r[i] = s[i];
				}
			}
			resid = Norm(r);
			if (print)
			{
				std::cout << "BICGSTAB " << std::setw(4) << iters << " " << std::setw(12) << resid << " | " << ftol;
				if( ptrue )
					std::cout << " true " << std::setw(12) << Resid(A, b, x);
				//std::cout << "\r";
				std::cout << std::endl;
				std::cout.flush();
			}
			iters++;
		} while ( resid > ftol && resid < dtol && omega && iters < maxiters+1);
		if (print)
		{
			std::cout << "BICGSTAB " << std::setw(4) << iters - 1 << " " << std::setw(12) << resid << " | " << ftol;
			if( ptrue )
				std::cout << " true " << Resid(A, b, x);
			//std::cout << "\r";
			std::cout << std::endl;
			std::cout.flush();
		}
		return ftol ? resid <= ftol : true;
	}
	size_t Bytes() const {return P.Bytes() + sizeof(const CSRMatrix *) + get_bytes(r) + get_bytes(r0) + get_bytes(y) + get_bytes(z) + get_bytes(p) + get_bytes(v) + get_bytes(s) + get_bytes(t);}
};

#endif //_BICGSTAB_H
