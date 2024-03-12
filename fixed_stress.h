#ifndef _FIXED_STRESS_H
#define _FIXED_STRESS_H
#include "method.h"
#include "vector.h"
#include "bicgstab.h"
#include "two_stage_gauss_seidel.h"
/*
 * Fixed-stress
 * this is fixed-stress method.
 * 
 * 
 * Implemented following
 * [1] A two-stage preconditioner for multiphase poromechanics in reservoir simulation
 * by White et al
 * 
 * 
 * 
 * Let the system (A x = b) be
 * 
 * | App   Apu | | p |   | bp |
 * |           | |   | = |    |
 * | Aup   Auu | | u |   | bu |
 * 
 * We compute the App shift by Dpp = -diag(Apu Auu^{-1} Aup e)
 * where e is a unit vector e = [1,...,1]^T.
 * 
 * The shifted and decoupled system is
 * | App+Dpp    0  | | p |   | bp |
 * |               | |   | = |    |
 * | Aup       Auu | | u |   | bu |
 * 
 * Solved by block Guass-Seidel method:
 * 1. p = (App + Dpp)^{-1} bpp
 * 2. u = (Aup)^{-1} (bu - Aup * p) 
 * 
 */

 
template<typename PressureSolver, typename DisplacementSolver>
class FixedStress : public virtual Methods, public TwoStageGaussSeidel<PressureSolver, DisplacementSolver>
{
	CSRMatrix B;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret = TwoStageGaussSeidel<PressureSolver, DisplacementSolver>::DefaultParameters();
		ret.Set("name","FixedStress");
		ret.Set("block_beg",0);
		ret.Set("block_end",0);
		ret.Set("type","fixed-stress # variants are fixed-stress, fixed-strain");
		ret.Set("check",0);
		ret.Set("order",0);
		ret.Set("verbosity",0);
		ret.Set("write_matrix",0);
		return ret;
	}
	FixedStress() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		std::string type = GetParameters().template Get<std::string>("type");
		bool print = GetParameters().template Get<int>("verbosity") ? true : false;
		bool check = GetParameters().template Get<int>("check") ? true : false;
		bool order = GetParameters().template Get<int>("order") ? true : false;
		bool write_matrix = GetParameters().template Get<int>("write_matrix") ? true : false;
		idx_t block_beg = GetParameters().template Get<idx_t>("block_beg");
		idx_t block_end = GetParameters().template Get<idx_t>("block_end");
		if( print )
			std::cout << "Method " << type << " pressure block " << block_beg << ":" << block_end << " write " << (write_matrix ? "yes" : "no") << std::endl;
		if( type == "fixed-stress" )
		{
			idx_t block_size = block_end - block_beg;
			std::vector<double> Duu(A.Size() - block_size, 0.0), Dpp(block_size, 0.0); //diagonal
			//extract inverse diagonal of the second block (displacement, Auu)
			{
				for (idx_t i = 0; i < block_beg; ++i)
				{
					for (idx_t j = 0; j < A.RowSize(i); ++j) 
						if (A.Col(i, j) == i) Duu[i] = A.Val(i, j);
				}
				for (idx_t i = block_end; i < A.Size(); ++i)
				{
					for (idx_t j = 0; j < A.RowSize(i); ++j)
						if (A.Col(i, j) == i) Duu[i - block_size] = A.Val(i, j);
				}
			}
			//multiplication Apu to Aup normalized by Duu
			//Go over rows of Apu
			for (idx_t i = block_beg; i < block_end; ++i)
			{
				for (idx_t j = 0; j < A.RowSize(i); ++j)
				{
					idx_t p = A.Col(i, j);
					double aip = A.Val(i, j);
					if (p < block_beg || p >= block_end)
					{
						//find element, contributing to diagonal in Aup
						for (idx_t l = 0; l < A.RowSize(p); ++l)
						{
							idx_t m = A.Col(p, l);
							double api = A.Val(p, l);
							if (m == i)
							{
								double duu = Duu[p >= block_end ? p - block_size : p];
								//std::cout << i << " add (" << i << "," << p << "," << aip << ") * (" << p << ',' << m << "," << api << ") / (" << (p >= block_end ? p - block_size : p) << "," << duu << ") = "  << aip * api / duu << std::endl;
								Dpp[i - block_beg] += aip * api / duu;
								//std::cout << "result: " << Dpp[i - block_beg] <<  " i " << i << " " << i - block_beg << std::endl;
							}
						}
					}
				}
			}
			//assemble shifted and decoupled matrix
			for (idx_t i = 0; i < block_beg; ++i)
			{
				for (idx_t j = 0; j < A.RowSize(i); ++j)
					B.PushBack(A.Col(i, j), A.Val(i, j));
				B.FinalizeRow();
			}
			for (idx_t i = block_beg; i < block_end; ++i)
			{
				for (idx_t j = 0; j < A.RowSize(i); ++j)
				{
					idx_t p = A.Col(i, j);
					if (p == i) //shift diagonal
						B.PushBack(p, A.Val(i, j) -Dpp[i - block_beg]);
					else if (p >= block_beg && p < block_end) //TODO!!!
						B.PushBack(p, A.Val(i, j));
					// zero off-diagonal block!
				}
				B.FinalizeRow();
			}
			for (idx_t i = block_end; i < A.Size(); ++i)
			{
				for (idx_t j = 0; j < A.RowSize(i); ++j)
					B.PushBack(A.Col(i, j), A.Val(i, j));
				B.FinalizeRow();
			}
			if (print)
			{
				std::cout << "Consumed by B: " << B.Bytes() / 1024.0 / 1024.0 << "Mb." << std::endl;
			}
			return TwoStageGaussSeidel<PressureSolver, DisplacementSolver>::Setup(B);
		}
		else if (type == "fixed-strain")
			return TwoStageGaussSeidel<PressureSolver,DisplacementSolver>::Setup(A);
		else
		{
			std::cout << "Error: unknown method type " << type << " expected either fixed-stress or fixed-strain." << std::endl;
			return false;
		}
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		return TwoStageGaussSeidel<PressureSolver,DisplacementSolver>::Solve(b,x);
	}
	size_t Bytes() const 
	{
		return B.Bytes() + TwoStageGaussSeidel<PressureSolver, DisplacementSolver>::Bytes();
	}
};


#endif //_FIXED_STRESS_H
