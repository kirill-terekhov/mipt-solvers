#ifndef _TWO_STAGE_GAUSS_SEIDEL_H
#define _TWO_STAGE_GAUSS_SEIDEL_H
#include "method.h"
/*
 * TwoStageGaussSeidel
 * This is the symmetric block gauss-seidel method for two blocks.
 */
template<typename Block1Solver, typename Block2Solver>
class TwoStageGaussSeidel : public virtual Methods
{
	const CSRMatrix * ptr_A;
	CSRMatrix B1, B2;
	Block1Solver iB1;
	Block2Solver iB2;
	mutable std::vector<double> b1,b2,x1,x2;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","TwoStageGaussSeidel");
		ret.Set("block_beg",0);
		ret.Set("block_end",0);
		ret.Set("check",0);
		ret.Set("verbosity",1);
		ret.Set("write_matrix",0);
		ret.SubParameters("BlockSolver") = Block1Solver::DefaultParameters();
		ret.SubParameters("SecondSolver") = Block2Solver::DefaultParameters();
		return ret;
	}
	TwoStageGaussSeidel() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		idx_t err = 0;
		ptr_A = &A;
		bool print = GetParameters().template Get<int>("verbosity") ? true : false;
		bool check = GetParameters().template Get<int>("check") ? true : false;
		bool write_matrix = GetParameters().template Get<int>("write_matrix") ? true : false;
		idx_t block_beg = GetParameters().template Get<idx_t>("block_beg");
		idx_t block_end = GetParameters().template Get<idx_t>("block_end");
		idx_t block_size = block_end-block_beg;
		if( print )
			std::cout << "System size " << A.Size() << " sub-system block " << block_beg << ":" << block_end << " size " << block_size << std::endl;
		if( block_size > A.Size() ) throw -1;
		
		//extract first block
		{
			for(idx_t i = block_beg; i < block_end; ++i)
			{
				for(idx_t j = 0; j < A.RowSize(i); ++j)
				{
					if( A.Col(i,j) >= block_beg && A.Col(i,j) < block_end )
						B1.PushBack(A.Col(i,j) - block_beg, A.Val(i,j));
				}
				B1.FinalizeRow();
			}
		}
		//extract second block
		{
			for(idx_t i = 0; i < block_beg; ++i)
			{
				for(idx_t j = 0; j < A.RowSize(i); ++j)
				{
					if( A.Col(i,j) < block_beg)
						B2.PushBack(A.Col(i,j), A.Val(i,j));
				}
				B2.FinalizeRow();
			}
			for(idx_t i = block_end; i < A.Size(); ++i)
			{
				for(idx_t j = 0; j < A.RowSize(i); ++j)
				{
					if( A.Col(i,j) >= block_end)
						B2.PushBack(A.Col(i,j)-block_size, A.Val(i,j));
				}
				B2.FinalizeRow();
			}
		}
		
		
		if( print )
		{
			std::cout << "Consumed by B1: " << B1.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by B2: " << B2.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
		}
		
		
		if( write_matrix )
		{
			if( print ) std::cout << "save B1.mtx" << std::endl;
			B1.Save("B1.mtx");
			if( print ) std::cout << "save B2.mtx" << std::endl;
			B2.Save("B2.mtx");
		}
		
		iB1.GetParameters() = GetParameters().SubParameters("BlockSolver");
		if( !iB1.Setup(B1) ) err++;
		iB2.GetParameters() = GetParameters().SubParameters("SecondSolver");
		if( !iB2.Setup(B2) ) err++;
		if( err )
		{
			std::cout << "Failed to setup block solver on " << err << " blocks." << std::endl;
			return false;
		}
		
		b1.resize(block_size);
		x1.resize(block_size);
		b2.resize(A.Size() - block_size);
		x2.resize(A.Size() - block_size);
				
		if( print )
		{
			std::cout << "Consumed by B1: " << B1.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by B2: " << B2.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by iB1: " << iB1.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by iB2: " << iB2.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by b1: " << get_bytes(b1)/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by x1: " << get_bytes(x1)/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by b2: " << get_bytes(b2)/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by x2: " << get_bytes(x2)/1024.0/1024.0 << "Mb." << std::endl;
		}
		return true;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool success = true;
		bool print = GetParameters().template Get<int>("verbosity") ? true : false;
		idx_t block_beg = GetParameters().template Get<idx_t>("block_beg");
		idx_t block_end = GetParameters().template Get<idx_t>("block_end");
		idx_t block_size = block_end - block_beg;
		const CSRMatrix & A = *ptr_A;
		x.resize(A.Size(),0.0);
		for(idx_t i = block_beg; i < block_end; ++i)
			x1[i-block_beg] = x[i];
		for(idx_t i = 0; i < block_beg; ++i)
		{
			b2[i] = b[i];
			x2[i] = x[i];
			//subtract updated off-second block part from residual
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p >= block_beg && p < block_end )
					b2[i] -= A.Val(i,j)*x1[p-block_beg];
			}
		}
		for(idx_t i = block_end; i < A.Size(); ++i)
		{
			b2[i-block_size] = b[i];
			x2[i-block_size] = x[i];
			//subtract updated off-second block part from residual
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p >= block_beg && p < block_end )
					b2[i-block_size] -= A.Val(i,j)*x1[p-block_beg];
			}
		}
		if( print ) std::cout << "Solve with the second block x2 " << Norm(x2) << " b2 " << Norm(b2) << std::endl;
		success &= iB2.Solve(b2, x2);
		for(idx_t i = block_beg; i < block_end; ++i)
		{
			b1[i-block_beg] = b[i];
			x1[i-block_beg] = x[i];
			//subtract updated off-first block part from residual
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p < block_beg )
					b1[i-block_beg] -= A.Val(i,j)*x2[p];
				else if( p >= block_end )
					b1[i-block_beg] -= A.Val(i,j)*x2[p-block_size];
			}
		}
		if( print ) std::cout << "Solve with the first block x1 " << Norm(x1) << " b1 " << Norm(b1) << std::endl;
		success &= iB1.Solve(b1, x1);
		if( print ) std::cout << " resid " << Resid(B1,b1,x1) << std::endl;
		if( !success ) 
		{
			std::cout << "Block 1 solver failed!" << std::endl;
			return false;
		}
		for(idx_t i = 0; i < block_beg; ++i)
		{
			b2[i] = b[i];
			x2[i] = x[i];
			//subtract updated off-second block part from residual
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p >= block_beg && p < block_end )
					b2[i] -= A.Val(i,j)*x1[p-block_beg];
			}
		}
		for(idx_t i = block_end; i < A.Size(); ++i)
		{
			b2[i-block_size] = b[i];
			x2[i-block_size] = x[i];
			//subtract updated off-second block part from residual
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p >= block_beg && p < block_end )
					b2[i-block_size] -= A.Val(i,j)*x1[p-block_beg];
			}
		}
		if( print ) std::cout << "solve with the second block x2 " << Norm(x2) << " b2 " << Norm(b2) << std::endl;
		success &= iB2.Solve(b2, x2);
		if( print ) std::cout << " resid " << Resid(B2,b2,x2) << std::endl;
		if( !success )
		{
			std::cout << "Block 2 solver failed!" << std::endl;
			return false;
		}
		//copy solution
		for(idx_t i = 0; i < block_beg; ++i)
			x[i] = x2[i];
		for(idx_t i = block_beg; i < block_end; ++i)
			x[i] = x1[i-block_beg];
		for(idx_t i = block_end; i < A.Size(); ++i)
			x[i] = x2[i-block_size];
		return success;
	}
	size_t Bytes() const 
	{
		size_t bytes = sizeof(CSRMatrix *) + B1.Bytes() + B2.Bytes();
		bytes += iB1.Bytes() + iB2.Bytes();
		bytes += get_bytes(b1) + get_bytes(x1) + get_bytes(b2) + get_bytes(x2);
		return bytes;
	}
};

#endif //_TWO_STAGE_GAUSS_SEIDEL_H
