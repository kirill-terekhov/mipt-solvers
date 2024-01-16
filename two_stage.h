#ifndef _TWO_STAGE_H
#define _TWO_STAGE_H
#include "method.h"
/*
 * TwoStage
 * This is the method that applies two preconditioners to the system.
 * The first preconditioner is applied to a (usually pressure) sub-block of the system.
 * The second preconditioner is applied to the full system with the account of the solution
 * of the first system in the right hand side.
 */
template<typename BlockSolver, typename SystemSolver>
class TwoStage : public virtual Methods
{
	const CSRMatrix * ptr_A;
	CSRMatrix B;
	BlockSolver iB;
	SystemSolver iA;
	mutable std::vector<double> bsub, xsub, bnew;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","TwoStage");
		ret.Set("block_beg",0);
		ret.Set("block_end",0);
		ret.Set("check",0);
		ret.Set("verbosity",1);
		ret.Set("write_matrix",0);
		ret.SubParameters("BlockSolver") = BlockSolver::DefaultParameters();
		ret.SubParameters("SecondSolver") = SystemSolver::DefaultParameters();
		return ret;
	}
	TwoStage() {GetParameters() = DefaultParameters();}
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
		
		
		//extract pressure block
		{
			for(idx_t i = block_beg; i < block_end; ++i)
			{
				for(idx_t j = 0; j < A.RowSize(i); ++j)
				{
					if( A.Col(i,j) >= block_beg && A.Col(i,j) < block_end )
						B.PushBack(A.Col(i,j) - block_beg, A.Val(i,j));
				}
				B.FinalizeRow();
			}
		}
		
		
		if( write_matrix )
		{
			if( print ) std::cout << "save B.mtx" << std::endl;
			B.Save("B.mtx");
		}
		
		iB.GetParameters() = GetParameters().SubParameters("BlockSolver");
		if( !iB.Setup(B) ) err++;
		iA.GetParameters() = GetParameters().SubParameters("SecondSolver");
		if( !iA.Setup(A) ) err++;
		if( err )
		{
			std::cout << "Failed to setup block solver on " << err << " blocks." << std::endl;
			return false;
		}
		
		bsub.resize(block_size);
		xsub.resize(block_size);
		bnew.resize(A.Size());
		
		if( print )
		{
			std::cout << "Consumed by B: " << B.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by iB: " << iB.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by iA: " << iA.Bytes()/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by bsub: " << get_bytes(bsub)/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by xsub: " << get_bytes(xsub)/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by bnew: " << get_bytes(bnew)/1024.0/1024.0 << "Mb." << std::endl;
		}
		return true;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool success = true;
		bool print = GetParameters().template Get<int>("verbosity") ? true : false;
		idx_t block_beg = GetParameters().template Get<idx_t>("block_beg");
		idx_t block_end = GetParameters().template Get<idx_t>("block_end");
		const CSRMatrix & A = *ptr_A;
		x.resize(A.Size(),0.0);
		for(idx_t i = block_beg; i < block_end; ++i)
		{
			bsub[i - block_beg] = b[i];
			xsub[i - block_beg] = x[i];
		}
		if( print )	std::cout << "Sub-system solution, x norm " << Norm(xsub) << " b norm " << Norm(bsub) << " residual " << Resid(B,bsub,xsub) << std::endl;
		success &= iB.Solve(bsub, xsub);
		if( print )	std::cout << "Residual " << Resid(B,bsub,xsub) << " solution norm " << Norm(xsub) << std::endl;
		if( !success ) 
		{
			std::cout << "Sub-system solver failed!" << std::endl;
			return false;
		}
		// A*(dx+xsub) = b
		// A*dx = b - A*xsub
		for(idx_t i = 0; i < A.Size(); ++i)
		{
			bnew[i] = b[i];
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p >= block_beg && p < block_end )
					bnew[i] -= A.Val(i,j)*xsub[p - block_beg];
			}
		}
		if( print )	std::cout << "Full-system solution, x norm " << Norm(x) << " modified b norm " << Norm(bnew) << " initial b norm " << Norm(b) << " residual " << Resid(A,bnew,x) << std::endl;
		success &= iA.Solve(bnew,x);
		if( print )	std::cout << "Residual " << Resid(A,bnew,x) << " solution norm " << Norm(x) << std::endl;
		if( !success )
		{
			std::cout << "Full system solver failed!" << std::endl;
			return false;
		}
		for(idx_t i = block_beg; i < block_end; ++i)
			x[i] += xsub[i - block_beg];
		if( print )
			std::cout << "Final residual " << Resid(A,b,x) << " solution norm " << Norm(x) << std::endl;
		return success;
	}
	size_t Bytes() const 
	{
		size_t bytes = sizeof(CSRMatrix *) + B.Bytes();
		bytes += iB.Bytes() + iA.Bytes();
		bytes += get_bytes(bsub) + get_bytes(xsub) + get_bytes(bnew);
		return bytes;
	}
};

#endif //_TWO_STAGE_H
