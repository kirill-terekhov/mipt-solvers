#ifndef _VANKA_H
#define _VANKA_H
#include "method.h"
#include "dense_matrix_routines.h"
#include <set>

#define SEQUENTIAL
/*
* Vanka smoother algorithm:
* collects overlapping blocks based on
* elements with zero diagonal
*/
#ifndef MAX_RAND
#define MAX_RAND 0x7fff
#endif
class Vanka : public Methods
{
	std::vector<double> D; //block diagonal
	std::vector<double> Dinv; //inverse of block diagonal
	std::vector<idx_t> bs; //block shifts in D and Dinv
	std::vector<idx_t> bia; //block sizes bia[i+1] - bia[i]
	std::vector<idx_t> bja; //block row indices
	mutable std::vector<double> r, g, y; //storage for intermediate result of current block
	const CSRMatrix * ptr_A;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","Vanka");
		ret.Set("omega",1.0);
		ret.Set("maxiters",2);
		ret.Set("verbosity", 1);
		ret.Set("check", 0);
		ret.Set("tol", 0);
		return ret;
	}
	Vanka() {GetParameters() = DefaultParameters();}
	void GetBlock(const CSRMatrix& A, const std::vector<idx_t>& rows, std::vector<double>& D)
	{
		D.resize(rows.size() * rows.size());
		std::fill(D.begin(), D.end(), 0.0);
		for (size_t l = 0; l < rows.size(); ++l)
			for (idx_t jt = 0; jt < A.RowSize(rows[l]); ++jt)
			{
				idx_t j = A.Col(rows[l], jt);
				for (idx_t q = 0; q < rows.size(); ++q) if (j == rows[q])
					D[l * rows.size() + q] = A.Val(rows[l], jt);
			}
	}
	void UnitMatrix(std::vector<double>& U, size_t s)
	{
		U.resize(s * s);
		std::fill(U.begin(), U.end(), 0.0);
		for (size_t l = 0; l < s; ++l)
			U[l * s + l] = 1.0;
	}
	bool ComputeBlocks(const CSRMatrix& A)
	{
		bool print = GetParameters().Get<int>("verbosity") & 1 ? true : false;
		std::vector<double> B, B0, iB, I;
		std::vector<double> diag(A.Size());
		std::vector< std::pair<double,idx_t> > row;
		std::vector<idx_t> rows;
		std::vector<int> order;
		std::map<size_t, int> sperb;
		size_t tot = 0, smin = std::numeric_limits<size_t>::max(), smax = 0, one = 1;
		A.Diagonal(diag);
		bia.push_back(0);
		for (idx_t k = 0; k < A.Size(); ++k)
		{
			bool done = false, bad = false;
			//zero diagonal
			if( 1 + diag[k] == 1)
			{
				//collect row elements without diagonal
				rows.push_back(k);
				for (idx_t j = 0; j < A.RowSize(k); ++j)
					if (A.Col(k, j) != k)
						rows.push_back(A.Col(k, j));
				//form and invert block
				GetBlock(A, rows, B);
				UnitMatrix(I, rows.size());
				iB.resize(rows.size() * rows.size());
				order.resize(rows.size());
				B0 = B;
				int ret = solvenxnxm(&B[0], &iB[0], &I[0], (int)rows.size(), (int)rows.size(), &order[0]);
				if (ret == 0)
				{
					transposenxn(&iB[0], (int)rows.size());
					D.insert(D.end(), B0.begin(), B0.end());
					Dinv.insert(Dinv.end(), iB.begin(), iB.end());
					for (size_t q = 0; q < rows.size(); ++q)
						bja.push_back(rows[q]);
					bs.push_back((idx_t)tot);
					tot += rows.size() * rows.size();
					smax = std::max(rows.size(), smax);
					smin = std::min(rows.size(), smin);
					sperb[rows.size()]++;
					done = true;

					bia.push_back((idx_t)bja.size()); //close block
				}
				rows.clear();
				row.clear();
			}
			
		}
		
		if (print)
		{
			std::cout << "Block size: " << smin << ":" << smax << " mean " << sqrt(tot / (bia.size() - 1.0)) << std::endl;
			std::cout << "Blocks: ";
			for (std::map<size_t, int>::iterator it = sperb.begin(); it != sperb.end(); ++it)
				std::cout << it->first << ":" << it->second << " ";
			std::cout << std::endl;
		}
		r.resize(smax);
		g.resize(smax);
		y.resize(smax);
		return true;
	}
	bool Setup(const CSRMatrix & A)
	{
		bool print = GetParameters().Get<int>("verbosity") & 1 ? true : false;
		bool success = true;
		idx_t err = 0;
		ptr_A = &A;
		success &= ComputeBlocks(A);
		if (print)
			std::cout << "Vanka consumed: " << Bytes() / 1024.0 / 1024.0 << " MB " << std::endl;
		return success;
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
			bool rev = iters % 2 ? true : false;
			for(idx_t it = 1; it < bia.size(); ++it)
			{
				idx_t i = rev ? (idx_t)bia.size() - it : it; //alternate backward and forward subsitution
				idx_t s = bia[i] - bia[i - 1];
				const double* bD = &D[bs[i - 1]];
				const double* bDinv = &Dinv[bs[i - 1]];
				for (idx_t k = bia[i - 1], q = 0; k < bia[i]; ++k, ++q)
				{
					idx_t j = bja[k]; // current row
					r[q] = b[j];
					for (idx_t jt = 0; jt < A.RowSize(j); ++jt)
						r[q] -= A.Val(j, jt) * x[A.Col(j, jt)];
					for (idx_t l = 0; l < s; ++l)
						r[q] += bD[q * s + l] * x[bja[bia[i - 1] + l]];
				}
				for (idx_t k = 0; k < s; ++k)
				{
					g[k] = 0.0;
					for (idx_t l = 0; l < s; ++l)
						g[k] += bDinv[k * s + l] * r[l];
				}
				bool print = false;
				for (idx_t k = bia[i - 1], q = 0; k < bia[i]; ++k, ++q)
				{
					idx_t j = bja[k]; // current row
					y[q] = omega * g[q] + (1.0 - omega) * x[j];
				}				
				for (idx_t k = bia[i - 1], q = 0; k < bia[i]; ++k, ++q)
				{
					idx_t j = bja[k]; // current row
					x[j] = y[q];
				}
			}
			if( tol )
			{
				resid = Resid(A,b,x);
				if( resid < tol ) break;
			}
			iters++;
		} while( iters < maxiters );
		if( print )
			std::cout << "Vanka " << iters << " true " << (tol ? resid : Resid(A,b,x)) << std::endl;
		return !tol || resid < tol;
	}
	size_t Bytes() const 
	{
		size_t bytes = sizeof(const CSRMatrix*);
		bytes += get_bytes(D) + get_bytes(Dinv) + get_bytes(bs);
		bytes += get_bytes(bia) + get_bytes(bja);
		bytes += get_bytes(r) + get_bytes(g) + get_bytes(y);
		return bytes;
	}
};

#endif //_VANKA_H
