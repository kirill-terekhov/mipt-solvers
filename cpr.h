#ifndef _CPR_H
#define _CPR_H
#include "method.h"
#include "vector.h"
/*
 * CPR
 * this is constrained pressure residual method.
 * 
 * Warning: the method detects the diagonal of the block, corresponding to saturations
 *          by negative value, which is true for water. However, not true for other phases.
 * 
 * TODO: detect sign by considering row-sum?
 * 
 * Implemented following
 * [1] A Constrained Pressure Residual Multiscale (CPR-MS) Compositional Solver
 * by Cusini et al
 * [2] Decoupling preconditioners in the implicit parallel accurate reservoir simulator (IPARS)
 * by Lacroix et al
 * 
 * 
 * 
 * Let the system (A x = b) be
 * 
 * | App   Aps | | p |   | bp |
 * |           | |   | = |    |
 * | Asp   Ass | | s |   | bs |
 * 
 * the left scaling coefficient is
 * 
 * |  I   -Dps Dss^{-1} |
 * |                    |
 * |            I       |
 * 
 * where Dps = diag(Aps) (a) = colsum(Aps) (b) = Aps (c)
 *       Dss = diag(Ass) (a) = colsum(Ass) (b)
 * 
 * since Dps may not be square, it is assumed that 
 * the first order upstream is used and the diagonal
 * terms are negative.
 * 
 * multiplying by the coefficient results in the transformation
 * 
 * | App - Dps Dss^{-1} Asp     Aps - Dps Dss^{-1} Ass | | p |   | bp - Dps Dss^{-1} bs |
 * |                                                   | |   | = |                      |
 * |         Asp                        Ass            | | s |   |          bs          |
 * 
 * where Aps - Dps Dss^{-1} Ass = 0 is assumed
 * 
 * The solution is split into two steps
 * 
 * 1. Solve pressure system
 * 
 * (App - Dps Dss^{-1} Asp) p1 = bp - Dps Dss^{-1} bs
 * 
 * 2. Solve full system
 * 
 * x = A^{-1} ( b - A [p1,0]^T ) + [p1,0]^T
 * 
 * 
 */
 
static CSRMatrix CPRScalingMatrix(const CSRMatrix & A, idx_t block_beg, idx_t block_end, bool colsum = true, double S_diag_sign = -1.0)
{
	CSRMatrix asm_C;
	idx_t block_size = (block_end-block_beg);
	std::vector<double> Q(A.Size() - block_size,0.0);
	if( colsum )
	{
		std::vector<double> D(A.Size() - block_size,0.0);
		for(idx_t i = block_beg; i < block_end; ++i)
		{
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p < block_beg )
					Q[p] += A.Val(i,j);
				else if( p >= block_end )
					Q[p-block_size] += A.Val(i,j);
			}
		}
		for(idx_t i = 0; i < block_beg; ++i)
		{
			for(idx_t j = 0; j < A.RowSize(i); ++j) 
			{
				idx_t p = A.Col(i,j);
				if( p < block_beg )
					D[p] += A.Val(i,j);
				else if( p >= block_end )
					D[p-block_size] += A.Val(i,j);
			}
		}
		for(idx_t i = block_end; i < A.Size(); ++i)
		{
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p < block_beg )
					D[p] += A.Val(i,j);
				else if( p >= block_end )
					D[p-block_size] += A.Val(i,j);
			}
		}	
		for(idx_t i = 0; i < A.Size() - block_size; ++i)
			Q[i] /= D[i];
	}
	else
	{
		for(idx_t i = block_beg; i < block_end; ++i)
		{
			for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j)*S_diag_sign > 0.0 ) //negative terms are diagonal
			{
				idx_t p = A.Col(i,j);
				if( p < block_beg ) 
					Q[p] = A.Val(i,j);
				else if( p >= block_end ) 
					Q[p-block_size] = A.Val(i,j);
			}
		}
		for(idx_t i = 0; i < block_beg; ++i)
		{
			for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Col(i,j) == i )
			{
				Q[i] /= A.Val(i,j);
				break;
			}
		}
		for(idx_t i = block_end; i < A.Size(); ++i)
		{
			for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Col(i,j) == i )
			{
				Q[i-block_size] /= A.Val(i,j);
				break;
			}
		}
	}
	for(idx_t i = 0; i < block_beg; ++i)
	{
		asm_C.PushBack(i,1.0);
		asm_C.FinalizeRow();
	}
	for(idx_t i = block_beg; i < block_end; ++i)
	{
		asm_C.PushBack(i,1.0); //unit diagonal
		//other pressure blocks, take out diagonal positions
		for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j)*S_diag_sign > 0.0 ) //negative terms outside of p-block are diagonal
		{
			idx_t p = A.Col(i,j);
			if( p < block_beg ) //check that we are outside of p-block
				asm_C.PushBack(p,-Q[p]);
			else if( p >= block_end )
				asm_C.PushBack(p,-Q[p-block_size]);
		}
		asm_C.FinalizeRow();
	}
	for(idx_t i = block_end; i < A.Size(); ++i)
	{
		asm_C.PushBack(i,1.0);
		asm_C.FinalizeRow();
	}
	if( false )
	{
		std::vector<int> ndiag(A.Size()-block_size,0);
		for(idx_t i = block_beg; i < block_end; ++i)
		{
			for(idx_t j = 0; j < asm_C.RowSize(i); ++j)
				if( asm_C.Col(i,j) < block_beg )
					ndiag[asm_C.Col(i,j)]++;
				else if( asm_C.Col(i,j) >= block_end )
					ndiag[asm_C.Col(i,j)-block_size]++;
		}
		int err = 0;
		for(idx_t i = 0; i < A.Size() - block_size; ++i)
			if( ndiag[i] > 1 )
			{
				//~ std::cout << ndiag[i] << " diagonals on line " << i << std::endl;
				err++;
			}
		if( err )
		{
			std::cout << "block: " << block_beg << ":" << block_end << std::endl;
			std::cout << "Save A.mtx" << std::endl;
			A.Save("A.mtx");
			std::cout << "Save C.mtx" << std::endl;
			asm_C.Save("C.mtx");
		}
	}
	return asm_C;
}
 
template<typename PressureSolver, typename SystemSolver, template<class,class> class TwoStageMethod>
class CPR : public virtual Methods, public TwoStageMethod<PressureSolver,SystemSolver>
{
	CSRMatrix C, CA;
	mutable std::vector<double> Cb;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret = TwoStageMethod<PressureSolver,SystemSolver>::DefaultParameters();
		ret.Set("name","CPR");
		ret.Set("block_beg",0);
		ret.Set("block_end",0);
		ret.Set("type","quasi # variants are none (no rescaling), quasi (quasi-IMPES), true (true IMPES)");
		ret.Set("check",0);
		ret.Set("order",0);
		ret.Set("verbosity",0);
		ret.Set("write_matrix",0);
		return ret;
	}
	CPR() {GetParameters() = DefaultParameters();}
	//colsum = true corresponds to True-IMPES, = false to Quasi-IMPES
	bool Setup(const CSRMatrix & A)
	{
		std::string type = GetParameters().Get<std::string>("type");
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		bool check = GetParameters().Get<int>("check") ? true : false;
		bool order = GetParameters().Get<int>("order") ? true : false;
		bool write_matrix = GetParameters().Get<int>("write_matrix") ? true : false;
		idx_t block_beg = GetParameters().Get<idx_t>("block_beg");
		idx_t block_end = GetParameters().Get<idx_t>("block_end");
		bool colsum = false;
		C.Clear();
		CA.Clear();
		if( print )
			std::cout << "CPR type " << type << " block " << block_beg << ":" << block_end << " write " << (write_matrix ? "yes" : "no") << std::endl;
		if( type != "none" )
		{
			if( type == "true" )
				colsum = true;
			else if( type == "quasi" )
				colsum = false;
			else
			{
				std::cout << "Error: unknown CPR type " << type << " expected either none or true or quasi." << std::endl;
				return false;
			}
			C = CPRScalingMatrix(A,block_beg,block_end,colsum);
			CA = C*A;
			if( order )
				CA.SortRows();
			Cb.resize(A.Size());
			if( write_matrix )
			{
				if( print ) std::cout << "save C.mtx" << std::endl;
				C.Save("C.mtx");
				if( print ) std::cout << "save CA.mtx" << std::endl;
				CA.Save("CA.mtx");
			}
			if( print )
			{
				std::cout << "Consumed by C: " << C.Bytes() / 1024.0 / 1024.0 << "Mb." << std::endl;
				std::cout << "Consumed by CA: " << CA.Bytes() / 1024.0 / 1024.0 << "Mb." << std::endl;
				std::cout << "Consumed by Cb: " << get_bytes(Cb) / 1024.0 / 1024.0 << "Mb." << std::endl;
			}
			return TwoStageMethod<PressureSolver,SystemSolver>::Setup(CA);
		}
		else return TwoStageMethod<PressureSolver,SystemSolver>::Setup(A);
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		bool scaling = GetParameters().Get<std::string>("type") == "none" ? false : true;
		if( scaling )
		{
			C.Multiply(1.0,b,0.0,Cb); //~ Cb = C*b
			return TwoStageMethod<PressureSolver,SystemSolver>::Solve(Cb,x);
		}
		else return TwoStageMethod<PressureSolver,SystemSolver>::Solve(b,x);
	}
	size_t Bytes() const 
	{
		return C.Bytes() + CA.Bytes() + get_bytes(Cb) + TwoStageMethod<PressureSolver,SystemSolver>::Bytes();
	}
};


template<typename Solver>
class CPRScaling : public virtual Methods
{
	CSRMatrix C, CA;
	mutable std::vector<double> Cb;
	Solver S;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","CPR");
		ret.Set("block_beg",0);
		ret.Set("block_end",0);
		ret.Set("type","quasi # variants are none (no rescaling), quasi (quasi-IMPES), true (true IMPES)");
		ret.Set("check",0);
		ret.Set("order",0);
		ret.Set("verbosity",0);
		ret.Set("write_matrix",0);
		ret.SubParameters("Solver") = Solver::DefaultParameters();
		return ret;
	}
	CPRScaling() {GetParameters() = DefaultParameters();}
	//colsum = true corresponds to True-IMPES, = false to Quasi-IMPES
	bool Setup(const CSRMatrix & A)
	{
		std::string type = GetParameters().Get<std::string>("type");
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		bool check = GetParameters().Get<int>("check") ? true : false;
		bool order = GetParameters().Get<int>("order") ? true : false;
		bool write_matrix = GetParameters().Get<int>("write_matrix") ? true : false;
		idx_t block_beg = GetParameters().Get<idx_t>("block_beg");
		idx_t block_end = GetParameters().Get<idx_t>("block_end");
		bool colsum = false;
		C.Clear();
		CA.Clear();
		if( print )
			std::cout << "CPR type " << type << " block " << block_beg << ":" << block_end << " write " << (write_matrix ? "yes" : "no") << std::endl;
		S.SetParameters(GetParameters().SubParameters("Solver"));
		if( type != "none" )
		{
			if( type == "true" )
				colsum = true;
			else if( type == "quasi" )
				colsum = false;
			else
			{
				std::cout << "Error: unknown CPR type " << type << " expected either none or true or quasi." << std::endl;
				return false;
			}
			C = CPRScalingMatrix(A,block_beg,block_end,colsum);
			CA = C*A;
			if( order )
				CA.SortRows();
			if( write_matrix )
			{
				if( print ) std::cout << "save C.mtx" << std::endl;
				C.Save("C.mtx");
				if( print ) std::cout << "save CA.mtx" << std::endl;
				CA.Save("CA.mtx");
			}
			bool success = S.Setup(CA);
			Cb.resize(A.Size());
			if( print )
			{
				std::cout << "Consumed by C: " << C.Bytes() / 1024.0 / 1024.0 << "Mb." << std::endl;
				std::cout << "Consumed by CA: " << CA.Bytes() / 1024.0 / 1024.0 << "Mb." << std::endl;
				std::cout << "Consumed by Cb: " << get_bytes(Cb) / 1024.0 / 1024.0 << "Mb." << std::endl;
			}
			return success;
		}
		else return S.Setup(A);
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		bool scaling = GetParameters().Get<std::string>("type") == "none" ? false : true;
		if( scaling )
		{
			C.Multiply(1.0,b,0.0,Cb); //~ Cb = C*b
			return S.Solve(Cb,x);
		}
		else return S.Solve(b,x);
	}
	size_t Bytes() const 
	{
		return C.Bytes() + CA.Bytes() + get_bytes(Cb) + S.Bytes();
	}
};


template<typename PressureSolver, typename SystemSolver>
class CPR2 : public Methods
{
	const CSRMatrix * ptr_A;
	std::vector<double> Dps;
	CSRMatrix Bpp;
	PressureSolver iBpp;
	SystemSolver iA;
	mutable std::vector<double> bp, bs, p1, r;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","CPR2");
		ret.Set("block_beg",0);
		ret.Set("block_end",0);
		ret.Set("type","quasi # variants are none (no rescaling), quasi (quasi-IMPES), true (true IMPES)");
		ret.Set("check",0);
		ret.Set("verbosity",0);
		ret.Set("write_matrix",0);
		ret.SubParameters("BlockSolver") = PressureSolver::DefaultParameters();
		ret.SubParameters("SecondSolver") = SystemSolver::DefaultParameters();
		return ret;
	}
	CPR2() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		CSRMatrix Asp;
		idx_t err = 0;
		ptr_A = &A;
		std::string type = GetParameters().Get<std::string>("type");
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		bool check = GetParameters().Get<int>("check") ? true : false;
		bool write_matrix = GetParameters().Get<int>("write_matrix") ? true : false;
		idx_t block_beg = GetParameters().Get<idx_t>("block_beg");
		idx_t block_end = GetParameters().Get<idx_t>("block_end");
		if( print )
			std::cout << "Matrix size " << A.Size() << " pressure block start " << block_beg << " end " << block_end << " size " << block_end-block_beg << " cpr type: " << type <<  std::endl;
		
		Bpp.Clear();
		Dps.clear();
		
		//assemble Dps
		if( type != "none" )
		{
			std::vector<double> Dss(A.Size() - (block_end - block_beg), 0.0); //diagonal or colsum
			Dps.resize(A.Size() - (block_end - block_beg),0.0);
			if( type == "true" )
			{
				for(idx_t i = block_beg; i < block_end; ++i)
				{
					for(idx_t j = 0; j < A.RowSize(i); ++j)
					{
						if( A.Col(i,j) < block_beg )
							Dps[A.Col(i,j)] += A.Val(i,j);
						else if( A.Col(i,j) >= block_end )
							Dps[A.Col(i,j)-block_end+block_beg] += A.Val(i,j);
					}
				}
				for(idx_t i = 0; i < block_beg; ++i)
				{
					for(idx_t j = 0; j < A.RowSize(i); ++j) 
					{
						if( A.Col(i,j) < block_beg )
							Dss[A.Col(i,j)] += A.Val(i,j);
						else if( A.Col(i,j) >= block_end )
							Dss[A.Col(i,j)-block_end+block_beg] += A.Val(i,j);
					}
				}
				for(idx_t i = block_end; i < A.Size(); ++i)
				{
					for(idx_t j = 0; j < A.RowSize(i); ++j)
					{
						if( A.Col(i,j) < block_beg )
							Dss[A.Col(i,j)] += A.Val(i,j);
						else if( A.Col(i,j) >= block_end )
							Dss[A.Col(i,j)-block_end+block_beg] += A.Val(i,j);
					}
				}	
			}
			else if( type == "quasi" )
			{
				for(idx_t i = block_beg; i < block_end; ++i)
				{
					for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j) < 0.0 ) //negative terms are diagonal
					{
						if( A.Col(i,j) < block_beg ) 
							Dps[A.Col(i,j)] = A.Val(i,j);
						else if( A.Col(i,j) >= block_end ) 
							Dps[A.Col(i,j)-block_end+block_beg] = A.Val(i,j);
					}
				}
				for(idx_t i = 0; i < block_beg; ++i)
				{
					for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Col(i,j) == i )
					{
						Dss[i] = A.Val(i,j);
						break;
					}
				}
				for(idx_t i = block_end; i < A.Size(); ++i)
				{
					for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Col(i,j) == i )
					{
						Dss[i-block_end+block_beg] = A.Val(i,j);
						break;
					}
				}
			}
			else
			{
				std::cout << "Unknown CPR type " << type << " expected \"quasi\" or \"true\"." << std::endl;
				return false;
			}
			
			if( print )
			{
				std::cout << "Consumed by Dps: " << get_bytes(Dps)/1024.0/1024.0 << "Mb." << std::endl;
				std::cout << "Consumed by Dss: " << get_bytes(Dss)/1024.0/1024.0 << "Mb." << std::endl;
			}
			if( write_matrix )
			{
				if( print ) std::cout << "save Dps.mtx" << std::endl;
				SaveVector("Dps.mtx",Dps);
				if( print ) std::cout << "save Dss.mtx" << std::endl;
				SaveVector("Dss.mtx",Dss);
			}
		
			if( check )
			{
				for(idx_t i = 0; i < Dss.size(); ++i) 
					if( !Dss[i] ) err++;
				if( err )
				{
					std::cout << "Error: " << err << " zeroes in Dss" << std::endl;
					return false;
				}
			}
				
			//scale Dps by Dss^{-1}
			for(idx_t i = 0; i < Dss.size(); ++i) 
				Dps[i] /= Dss[i];
			
			if( write_matrix )
			{
				if( print ) std::cout << "save DpsDss.mtx" << std::endl;
				SaveVector("DpsDss.mtx",Dps);
			}
		}
		//assemble Bpp
		{
			for(idx_t i = block_beg; i < block_end; ++i)
			{
				for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j) )
				{
					if( A.Col(i,j) >= block_beg && A.Col(i,j) < block_end )
						Bpp.PushBack(A.Col(i,j) - block_beg, A.Val(i,j));
				}
				Bpp.FinalizeRow();
			}
		}
		//assemble Asp
		{
			for(idx_t i = 0; i < block_beg; ++i)
			{
				for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j) )
				{
					if( A.Col(i,j) >= block_beg && A.Col(i,j) < block_end )
						Asp.PushBack(A.Col(i,j)-block_beg,A.Val(i,j));
				}
				Asp.FinalizeRow();
			}
			for(idx_t i = block_end; i < A.Size(); ++i)
			{
				for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j) )
				{
					if( A.Col(i,j) >= block_beg && A.Col(i,j) < block_end )
						Asp.PushBack(A.Col(i,j)-block_beg,A.Val(i,j));
				}
				Asp.FinalizeRow();
			}
		}
		
		if( print )
		{
			std::cout << "Consumed by App: " << Bpp.Bytes()/1024.0/1024.0 << "Mb nonzeros " << Bpp.Nonzeros() << "." << std::endl;
			std::cout << "Consumed by Asp: " << Asp.Bytes()/1024.0/1024.0 << "Mb nonzeros " << Asp.Nonzeros() << "." << std::endl;
		}
		
		if( check )
		{
			for(idx_t i = block_beg; i < block_end; ++i)
			{
				for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j) )
				{
					if( A.Col(i,j) < block_beg )
						Aps.PushBack(A.Col(i,j),A.Val(i,j));
					else if( A.Col(i,j) >= block_end )
						Aps.PushBack(A.Col(i,j)-block_end+block_beg,A.Val(i,j));
				}
				Aps.FinalizeRow();
			}
			CSRMatrix Aps = asm_Aps;
			
			std::cout << "A norm: " << Norm(A.get_a()) << std::endl;
			std::cout << "App norm: " << Norm(Bpp.get_a()) << std::endl;
			std::cout << "Aps norm: " << Norm(Aps.get_a()) << std::endl;
			std::cout << "Asp norm: " << Norm(Asp.get_a()) << std::endl;
				
			if( type != "none" )
			{
				CSRMatrix asm_Ass;
				for(idx_t i = 0; i < block_beg; ++i)
				{
					for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j) )
					{
						if( A.Col(i,j) < block_beg )
							asm_Ass.PushBack(A.Col(i,j),A.Val(i,j));
						else if( A.Col(i,j) >= block_end )
							asm_Ass.PushBack(A.Col(i,j)-block_end+block_beg,A.Val(i,j));
					}
					asm_Ass.FinalizeRow(i);
				}
				for(idx_t i = block_end; i < A.Size(); ++i)
				{
					asm_Ass.InitializeRow(i+block_beg-block_end);
					for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Val(i,j) )
					{
						if( A.Col(i,j) < block_beg )
							asm_Ass.PushBack(A.Col(i,j),A.Val(i,j));
						else if( A.Col(i,j) >= block_end )
							asm_Ass.PushBack(A.Col(i,j)-block_end+block_beg,A.Val(i,j));
					}
					asm_Ass.FinalizeRow(i+block_beg-block_end);
				}
				CSRMatrix Ass = asm_Ass;
				CSRMatrix Bps = Aps - Scale(Dps,Ass);
				std::cout << "Ass norm: " << Norm(Ass.get_a()) << std::endl;
				std::cout << "Bps norm: " << Norm(Bps.get_a()) << std::endl;
				if( write_matrix )
				{
					if( print ) std::cout << "save Aps.mtx" << std::endl;
					Aps.Save("Aps.mtx");
					if( print ) std::cout << "save ApsT.mtx" << std::endl;
					Aps.Transpose().Save("ApsT.mtx");
					if( print ) std::cout << "save Bps.mtx" << std::endl;
					Bps.Save("Bps.mtx");
					if( print ) std::cout << "save Ass.mtx" << std::endl;
					Ass.Save("Ass.mtx");
					if( print ) std::cout << "save AssT.mtx" << std::endl;
					Ass.Transpose().Save("AssT.mtx");
				}
			}
		}
		
		if( write_matrix )
		{
			if( print ) std::cout << "save App.mtx" << std::endl;
			Bpp.Save("App.mtx");
		}
		
		if( type != "none" )
			Bpp = Bpp - Scale(Dps,Asp);
		Bpp.SortRows();
			
		if( check )
			std::cout << "Bpp norm: " << Norm(Bpp.get_a()) << std::endl;
		
		if( print )
			std::cout << "Consumed by Bpp: " << Bpp.Bytes()/1024.0/1024.0 << "Mb nonzeros " << Bpp.Nonzeros() << "." << std::endl;
		
		if( write_matrix )
		{
			if( print ) std::cout << "save Bpp.mtx" << std::endl;
			Bpp.Save("Bpp.mtx");
			if( print ) std::cout << "save Asp.mtx" << std::endl;
			Asp.Save("Asp.mtx");
		}
		
		iBpp.GetParameters() = GetParameters().SubParameters("BlockSolver");
		if( !iBpp.Setup(Bpp) ) err++;
		iA.GetParameters()   = GetParameters().SubParameters("SecondSolver");
		if( !iA.Setup(A) ) err++;
		if( err )
		{
			std::cout << "Failed to setup block solver on " << err << " blocks." << std::endl;
			return false;
		}
		
		bp.resize(block_end - block_beg);
		bs.resize(A.Size() - (block_end - block_beg));
		p1.resize(block_end - block_beg);
		r.resize(A.Size());
				
		if( print )
		{
			std::cout << "Consumed by bp: " << bp.capacity()*sizeof(double)/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by bs: " << bs.capacity()*sizeof(double)/1024.0/1024.0 << "Mb." << std::endl;
			std::cout << "Consumed by p1: " << p1.capacity()*sizeof(double)/1024.0/1024.0 << "Mb." << std::endl;
		}
		
		//~ static int cnt = 0;
		//~ if( cnt++ == 1 ) exit(-1);
		return true;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		//~ return iA.Solve(b,x);
		bool success = true;
		bool print = GetParameters().Get<int>("verbosity") ? true : false;
		bool have_Dps = GetParameters().Get<std::string>("type") == "none" ? false : true;
		idx_t block_beg = GetParameters().Get<idx_t>("block_beg");
		idx_t block_end = GetParameters().Get<idx_t>("block_end");
		const CSRMatrix & A = *ptr_A;
		x.resize(A.Size(),0.0);
		{
			for(idx_t i = block_beg; i < block_end; ++i)
				bp[i - block_beg] = b[i];
			for(idx_t i = 0; i < block_beg; ++i)
				bs[i] = b[i];
			for(idx_t i = block_end; i < A.Size(); ++i)
				bs[i-block_end+block_beg] = b[i];
		}
		if( have_Dps )
		{
			//~ Dps.Multiply(-1.0,bs,1.0,bp); //~ bp -= Dps*bs
			for(idx_t i = block_beg; i < block_end; ++i)
				bp[i] -= Dps[i]*bs[i];
		}
		success &= iBpp.Solve(bp, p1);
		if( !success ) 
		{
			std::cout << "Pressure system solver failed!" << std::endl;
			return false;
		}
		for(idx_t i = 0; i < A.Size(); ++i)
		{
			r[i] = b[i];
			for(idx_t j = 0; j < A.RowSize(i); ++j)
			{
				idx_t p = A.Col(i,j);
				if( p >= block_beg && p < block_end )
					r[i] -= A.Val(i,j)*p1[p - block_beg];
			}
		}
		success &= iA.Solve(r,x);
		if( !success )
		{
			std::cout << "Full system solver failed!" << std::endl;
			return false;
		}
		for(idx_t i = block_beg; i < block_end; ++i)
			x[i] += p1[i - block_beg];
		return success;
	}
	size_t Bytes() const 
	{
		size_t bytes = sizeof(CSRMatrix *) + get_bytes(Dps) + Bpp.Bytes();
		bytes += iA.Bytes() + iBpp.Bytes();
		bytes += get_bytes(bp) + get_bytes(bs) + get_bytes(p1) + get_bytes(r);
		return bytes;
	}
};


#endif //_CPR_H
