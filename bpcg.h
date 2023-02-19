#ifndef _BPCG_H
#define _BPCG_H
#include "method.h"
#include <cmath>

/*
* Conjugate gradient method of Bramble-Pasciak
* [1] A Preconditioning Technique for Indefinite Systems Resulting from Mixed Approximations of Elliptic Problems
* James Bramble and Joseph Pasciak
*/

template<typename Preconditioner>
class BramblePasciakCG : public Methods
{
	Preconditioner P;
	mutable std::vector<double> r1, r2, p1, p2, q1, q2;
	mutable std::vector<double> z1, z2, f1, f2, x1, x2;
	mutable std::vector<double> w1, w2, pr, apr;
	std::vector<bool> excl;
	std::vector<int> off;
	const CSRMatrix * ptr_A;
	CSRMatrix B, E, F, C; //matrix splitting
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","BPCG");
		ret.Set("B_block_beg", 0);
		ret.Set("B_block_end", 0);
		ret.Set("tol",1.0e-10);
		ret.Set("rtol",1.0e-7);
		ret.Set("dtol",1.0e+10);
		ret.Set("maxiters",5000);
		ret.Set("verbosity",1);
		ret.Set("true_residual", 1);
		ret.Set("flip_sign", 0);
		ret.SubParameters("Preconditioner") = Preconditioner::DefaultParameters();
		return ret;
	}
	BramblePasciakCG() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		bool success;
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		int Bbeg = GetParameters().Get<int>("B_block_beg");
		int Bend = GetParameters().Get<int>("B_block_end");
		int flip = GetParameters().Get<int>("flip_sign");
		if (Bbeg == Bend)
		{
			std::cout << "No block selected!" << std::endl;
			return false;
		}
		P.SetParameters(GetParameters().SubParameters("Preconditioner"));
		ptr_A = &A;
		if( print ) std::cout << "A symmetric: " << (A.Symmetric() ? "yes" : "no") << " size " << A.Size() << " nnz " << A.Nonzeros() << std::endl;
		if(print) std::cout << "B block: " << Bbeg << ":" << Bend << std::endl;
		//first mark dirichlet conditions (they could break the symmetry)
		excl.resize(A.Size(), false);
		off.resize(A.Size(), 0);
		int nskip = 0;
		
		for (idx_t k = 0; k < A.Size(); ++k)
		{
			off[k] = nskip;
			if (A.RowSize(k) == 1 && A.Col(k, 0) == k)
			{
				excl[k] = true;
				nskip++;
			}
		}
		
		std::cout << "Skip dirichlet conditions: " << nskip << std::endl;
		//assemble blocks B F
		for (idx_t k = Bbeg; k < Bend; ++k) if (!excl[k])
		{
			for (idx_t l = 0; l < A.RowSize(k); ++l) if (!excl[A.Col(k, l)])
			{
				if (A.Col(k, l) >= Bbeg && A.Col(k, l) < Bend)
					B.PushBack(A.Col(k, l) - (off[A.Col(k, l)] - off[Bbeg]) - Bbeg, A.Val(k, l));
				else if (A.Col(k, l) < Bbeg)
					F.PushBack(A.Col(k, l) - off[A.Col(k, l)], A.Val(k, l));
				else
					F.PushBack(A.Col(k, l) - (off[A.Col(k, l)] - (off[Bend] - off[Bbeg])) - (Bend - Bbeg), A.Val(k, l));
			}
			B.FinalizeRow();
			F.FinalizeRow();
		}
		if (print)
		{
			std::cout << "B size " << B.Size() << " nnz " << B.Nonzeros() << std::endl;
			std::cout << "F size " << F.Size() << " nnz " << F.Nonzeros() << std::endl;
			std::cout << "B symmetric? " << (B.Symmetric() ? "yes" : "no") << " frobenius norm " << B.FrobeniusNorm() << std::endl;
		}
		//assemble blocks E C with minus sign
		for (idx_t k = 0; k < A.Size(); ++k) if ((k < Bbeg || k >= Bend) && !excl[k])
		{
			for (idx_t l = 0; l < A.RowSize(k); ++l) if (!excl[A.Col(k, l)])
			{
				if (A.Col(k, l) >= Bbeg && A.Col(k, l) < Bend)
					E.PushBack(A.Col(k, l) - (off[A.Col(k, l)] - off[Bbeg]) - Bbeg, (1 - 2 * flip) * A.Val(k, l));
				else if (A.Col(k, l) < Bbeg)
					C.PushBack(A.Col(k, l) - off[A.Col(k, l)], (1 - 2 * flip) * A.Val(k, l));
				else
					C.PushBack(A.Col(k, l) - (off[A.Col(k, l)] - (off[Bend] - off[Bbeg])) - (Bend - Bbeg), (1 - 2 * flip) * A.Val(k, l));
			}
			C.FinalizeRow();
			E.FinalizeRow();
		}
		if (print)
		{
			std::cout << "E size " << E.Size() << " nnz " << E.Nonzeros() << std::endl;
			std::cout << "C size " << C.Size() << " nnz " << C.Nonzeros() << std::endl;
			std::cout << "C symmetric? " << (C.Symmetric() ? "yes" : "no") << " frobenius norm " << C.FrobeniusNorm() << std::endl;
		}
		if (print)
		{
			std::cout << "||E - F^T||: " << (E - F.Transpose()).FrobeniusNorm() << std::endl;
			std::cout << "||E + F^T||: " << (E + F.Transpose()).FrobeniusNorm() << std::endl;
			E.Save("E.mtx");
			F.Save("F.mtx");
		}
		success = P.Setup(B);
		f1.resize(B.Size());
		f2.resize(E.Size());
		x1.resize(B.Size());
		x2.resize(E.Size());
		p1.resize(B.Size());
		p2.resize(E.Size());
		r1.resize(B.Size());
		r2.resize(E.Size());
		w1.resize(B.Size());
		w2.resize(E.Size());
		q1.resize(B.Size());
		q2.resize(E.Size());
		pr.resize(B.Size());
		apr.resize(B.Size());
		return success;
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
	void RecordSolution(int Bbeg, int Bend, const std::vector<double>& x1, const std::vector<double>& x2, std::vector<double>& x) const
	{
		for (idx_t k = Bbeg; k < Bend; ++k) if (!excl[k])
			x[k] = x1[k - (off[k] - off[Bbeg]) - Bbeg];
		for (idx_t k = 0; k < x.size(); ++k) if ((k < Bbeg || k >= Bend) && !excl[k])
		{
			if (k < Bbeg)
				x[k] = x2[k - off[k]];// *(1 - 2 * flip);
			else
				x[k] = x2[k - (off[k] - (off[Bend] - off[Bbeg])) - (Bend - Bbeg)];// *(1 - 2 * flip);
		}
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		int Bbeg = GetParameters().Get<int>("B_block_beg");
		int Bend = GetParameters().Get<int>("B_block_end");
		int flip = GetParameters().Get<int>("flip_sign");
		bool ptrue = GetParameters().Get<int>("true_residual") ? true : false;
		int maxiters = GetParameters().Get<int>("maxiters");
		double tol  = GetParameters().Get<double>("tol");
		double rtol = GetParameters().Get<double>("rtol");
		double dtol = GetParameters().Get<double>("dtol");
		int iters = 1;
		double resid, resid0, beta, alpha, kappa, wr, wrn, ftol;
		const CSRMatrix & A = *ptr_A;
		idx_t size = A.Size();
		x.resize(size, 0.0);
        //extract rhs with account of excluded entries
		//precompute Dirichlet values
		for (idx_t k = 0; k < A.Size(); ++k) if (excl[k]) 
			x[k] = b[k] / A.Val(k, 0);
		//eliminate dirichlet values from rhs
		//B and F part
		for (idx_t k = Bbeg; k < Bend; ++k) if (!excl[k])
		{
			x1[k - (off[k] - off[Bbeg]) - Bbeg] = x[k];
			f1[k - (off[k] - off[Bbeg]) - Bbeg] = b[k];
			for (idx_t l = 0; l < A.RowSize(k); ++l) if (excl[A.Col(k, l)])
				f1[k - (off[k] - off[Bbeg]) - Bbeg] -= A.Val(k, l) * x[A.Col(k, l)];
		}
		//E and C part
		for (idx_t k = 0; k < A.Size(); ++k) if ((k < Bbeg || k >= Bend) && !excl[k])
		{
			double vf2 = b[k];
			for (idx_t l = 0; l < A.RowSize(k); ++l) if (excl[A.Col(k, l)])
				vf2 -= A.Val(k, l) * x[A.Col(k, l)];
			if (k < Bbeg)
			{
				f2[k - off[k]] = vf2;
				x2[k - off[k]] = x[k];
			}
			else
			{
				f2[k - (off[k] - (off[Bend] - off[Bbeg])) - (Bend - Bbeg)] = (1 - 2 * flip) * vf2;
				x2[k - (off[k] - (off[Bend] - off[Bbeg])) - (Bend - Bbeg)] = x[k];
			}
		}
		// F = [ P f1 \\ E P f1 - f2 ]
		ApplyPreconditioner(f1, pr);
		std::copy(pr.begin(), pr.end(), p1.begin());
		std::copy(f2.begin(), f2.end(), p2.begin());
		E.Multiply(1.0, pr, -1.0, p2);
		
		std::copy(f1.begin(), f1.end(), r1.begin());
		std::copy(f2.begin(), f2.end(), r2.begin());
		B.Multiply(1.0, pr, -1.0, r1);
		E.Multiply(1.0, pr, -1.0, r2);
		
		kappa = Dot(r1, p1) + Dot(r2, p2);
		resid0 = resid = sqrt(fabs(kappa));
		ftol = std::max(tol,rtol*resid0);
		if (print)
		{
			std::cout << "BPCG " << std::setw(4) << 0 << " " << std::setw(12) << resid << " | " << ftol;
			if (ptrue)
			{
				RecordSolution(Bbeg,Bend,x1, x2, x);
				std::cout << " true " << std::setw(12) << Resid(A, b, x);
			}
			std::cout << std::endl;
		}
        while( resid > ftol && resid < dtol && iters < maxiters+1 )
		{
			B.Multiply(1.0, p1, 0.0, q1);
			F.Multiply(1.0, p2, 1.0, q1);
			ApplyPreconditioner(q1, apr);
			B.Multiply(1.0, apr, -1.0, q1);
			E.Multiply(1.0, apr, 0.0, q2);
			E.Multiply(-1.0, p1, 1.0, q2);

			alpha = kappa / (Dot(p1, q1) + Dot(p2, q2));

			for (idx_t k = 0; k < x1.size(); ++k) x1[k] += alpha * p1[k];
			for (idx_t k = 0; k < x2.size(); ++k) x2[k] += alpha * p2[k];

			for (idx_t k = 0; k < r1.size(); ++k) r1[k] -= alpha * q1[k];
			for (idx_t k = 0; k < r2.size(); ++k) r2[k] -= alpha * q2[k];

			for (idx_t k = 0; k < pr.size(); ++k) pr[k] -= alpha * apr[k];

			for (idx_t k = 0; k < w1.size(); ++k) w1[k] = pr[k];
			for (idx_t k = 0; k < w2.size(); ++k) w2[k] = r1[k];

			beta = 1.0 / kappa;
			kappa = Dot(w1, r1) + Dot(w2, r2);	
			beta *= kappa;

			for(idx_t k = 0; k < p1.size(); ++k) p1[k] = w1[k] + beta * p1[k];
			for(idx_t k = 0; k < p2.size(); ++k) p2[k] = w2[k] + beta * p2[k];

			resid = std::sqrt(fabs(kappa));
			if (print)
			{
				std::cout << "BPCG "  << std::setw(4) << iters << " " << std::setw(14) << resid << " | " << ftol;
				if (ptrue)
				{
					RecordSolution(Bbeg,Bend,x1, x2, x);
					std::cout << " true " << std::setw(14) << Resid(A, b, x);
				}
				std::cout << std::endl;
			}
			iters++;
		}
		RecordSolution(Bbeg,Bend,x1, x2, x);
		if (print)
		{
			std::cout << "BPCG " << std::setw(4) << iters;
			if( ptrue ) std::cout << " true " << std::setw(14) << Resid(A, b, x);
			std::cout << std::endl;
		}
		return resid <= ftol;
	}
	size_t Bytes() const 
	{
		size_t ret = sizeof(const CSRMatrix*);
		ret += P.Bytes();
		ret += B.Bytes();
		ret += F.Bytes();
		ret += E.Bytes();
		ret += get_bytes(f1) + get_bytes(f2);
		ret += get_bytes(x1) + get_bytes(x2);
		ret += get_bytes(r1) + get_bytes(r2);
		ret += get_bytes(z1) + get_bytes(z2);
		ret += get_bytes(p1) + get_bytes(p2);
		ret += get_bytes(q1) + get_bytes(q2);
		ret += get_bytes(w1) + get_bytes(w2);
		ret += get_bytes(pr) + get_bytes(apr);
		return ret;}
};

#endif //_BPCG_H
