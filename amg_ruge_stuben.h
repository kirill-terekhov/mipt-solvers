#ifndef _AMG_RUGE_STUBEN_H
#define _AMG_RUGE_STUBEN_H
#include "method.h"
#include "priority_queue.h"
/*
 * AMGRugeStuben
 * Implemented following
 * "Algebraic multigrid"
 * by J.W. Ruge and K. Stuben
 * 
 * There is a modification of algorithm A.3 on page 102
 * based on the following description in text.
 * 
 * The algorithm is split into two functions:
 * CFRefine
 * Operator
 */

template<typename Smoother, typename CoarsestSolver>
class AMGRugeStuben : public Methods
{
	idx_t CSize;
	Smoother G;
	CSRMatrix B, P, R;
	const CSRMatrix * ptr_A;
	Methods * Next;
	mutable std::vector<double> rn, xn, r;
	const char C, F, U;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","AMGRugeStuben");
		ret.Set("write_matrix",0);
		ret.Set("level","*");
		ret.Set("phi",0.25);
		ret.Set("verbosity",1);
		ret.Set("check",1);
		ret.Set("order",0);
		ret.Set("cycle","V"); // V, F, W
		ret.SubParameters("Smoother") = Smoother::DefaultParameters();
		ret.SubParameters("CoarsestSolver") = CoarsestSolver::DefaultParameters();
		ret.SubParameters("CoarsestSolver").Set("verbosity","0");
		return ret;
	}
	AMGRugeStuben() : Next(NULL), C(1), F(2), U(4) {GetParameters() = DefaultParameters();}
	~AMGRugeStuben() { if(Next) delete Next;}
	void StrongGraph(const CSRMatrix & A, const std::vector<double> & diag, double phi, CSRMatrixType<idx_t> & S)
	{
		S.Clear();
		for(idx_t i = 0; i < A.Size(); ++i)
		{
			double rowmax = 0, sign = -1;
			if( diag[i] < 0 ) sign = 1;
			for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Col(i,j) != i )
				rowmax = std::max(rowmax,sign*A.Val(i,j));
			for(idx_t j = 0; j < A.RowSize(i); ++j) if( A.Col(i,j) != i )
				if( sign*A.Val(i,j) > phi*rowmax ) S.PushBack(A.Col(i,j),j);
			S.FinalizeRow();
		}
	}
	void CFSplitting(const CSRMatrixType<idx_t> & S, const CSRMatrixType<idx_t> & St, std::vector<char> & CF) const
	{
		CF.resize(St.Size());
		std::fill(CF.begin(),CF.end(),U);
		PriorityQueue<idx_t,std::greater<idx_t> > queue(St.Size());
		for(idx_t i = 0; i < St.Size(); ++i)
			queue.Push(i,St.RowSize(i));
		while( !queue.Empty() )
		{
			idx_t i = queue.Pop();
			if( CF[i] == U )
			{
				CF[i] = C;
				for(idx_t jt = 0; jt < St.RowSize(i); ++jt) 
				{
					idx_t j = St.Col(i,jt);
					if( CF[j] == U )
					{
						CF[j] = F;
						for(idx_t lt = 0; lt < S.RowSize(j); ++lt)
						{
							idx_t l = S.Col(j,lt);
							if( CF[l] == U ) 
								queue.IncreaseKey(l,queue.GetKey(l)+1);
						}
					}
				}
				for(idx_t jt = 0; jt < S.RowSize(i); ++jt) 
				{
					idx_t j = S.Col(i,jt);
					if( CF[j] == U ) 
						queue.DecreaseKey(j,queue.GetKey(j)-1);
				}
			}
		}
	}
	void CFRefine(const CSRMatrix & A, const CSRMatrixType<idx_t> & S, std::vector<char> & CF) const
	{
		std::vector<idx_t> Cnew;
		RowAccumulator<char> list(A.Size());
		for(idx_t i = 0; i < A.Size(); ++i) 
		{
			if( CF[i] == F )
			{
				//find out weak and strong skipped values
				for(idx_t jt = 0; jt < A.RowSize(i); ++jt)
				{
					idx_t j = A.Col(i,jt);
					if( j != i ) list[j] = CF[j];
				}
				//mark strong
				for(idx_t jt = 0; jt < S.RowSize(i); ++jt)
				{
					idx_t j = S.Col(i,jt);
					list[j] |= U;
				}
				//iterate over strong F-points
				for(idx_t j = list.Begin(); j != list.End(); j = list.Next(j))
				{
					if( list.Get(j) == (F|U) ) // strong skipped
					{
						//intersect S[j] and (C|U)
						bool found = false;
						for(idx_t kt = 0; kt < S.RowSize(j); kt++)
						{
							idx_t k = S.Col(j,kt);
							if( list.Contains(k) && list.Get(k) == (C|U) )
							{
								found = true;
								break;
							}
						}
						if( !found ) Cnew.push_back(j);
					}
				}
				if( Cnew.size() > 1 )
					CF[i] = C;
				else if( !Cnew.empty() )
					CF[Cnew.front()] = C;
				list.Clear();
				Cnew.clear();
			}
		}
	}
	void Operator(const CSRMatrix & A, const CSRMatrixType<idx_t> & S, const std::vector<char> & CF, CSRMatrix & I) const
	{
		bool   check        = GetParameters().Get<int>("check") ? true : false;
		I.Clear();
		std::vector<idx_t> id;
		id.resize(A.Size(),-1);
		idx_t cnt = 0;
		for(idx_t i = 0; i < A.Size(); ++i) 
			if( CF[i] == C ) id[i] = cnt++;
		{
			double sum, wgt, rowsum, rownorm;
			RowAccumulator< std::pair<char,idx_t> > list(A.Size());
			RowAccumulator<double> vals(A.Size());
			for(idx_t i = 0; i < A.Size(); ++i) 
			{
				if( CF[i] == F )
				{
					//find out weak and strong skipped values
					rowsum = 0;
					for(idx_t jt = 0; jt < A.RowSize(i); ++jt)
					{
						idx_t j = A.Col(i,jt);
						if( j != i ) list[j] = std::make_pair(CF[j],jt);
						else vals[i] = A.Val(i,jt); //diagonal
						if( check ) rowsum += A.Val(i,jt);
					}
					//mark strong
					for(idx_t jt = 0; jt < S.RowSize(i); ++jt)
					{
						idx_t j = S.Col(i,jt);
						list.Get(j).first |= U; //strong
						if( list.Get(j).first == (C|U) ) //interpolatory
							vals[j] = A.Val(i,list.Get(j).second);
					}
					for(idx_t j = list.Begin(); j != list.End(); j = list.Next(j))
					{
						if( list.Get(j).first == (F|U) ) //strong non-inpolatory
						{
							sum = 0;
							for(idx_t kt = 0; kt < S.RowSize(j); ++kt)
							{
								idx_t k = S.Col(j,kt);
								if( list.Contains(k) && list.Get(k).first == (C|U) ) //interpolatory
									sum += A.Val(j,S.Val(j,kt));
							}
							if( fabs(sum) > 1.0e-7 )
							{
								for(idx_t kt = 0; kt < S.RowSize(j); ++kt)
								{
									idx_t k = S.Col(j,kt);
									if( list.Contains(k) && list.Get(k).first == (C|U) ) //interpolatory
										vals.Get(k) += A.Val(i,list.Get(j).second)*A.Val(j,S.Val(j,kt))/sum;
								}
							}
							else vals.Get(i) += A.Val(i,list.Get(j).second);
						}
						else if( list.Get(j).first != (C|U) ) //weak non-interpolatory
							vals.Get(i) += A.Val(i,list.Get(j).second);
					}
					if( check )
					{
						rowsum -= vals.Get(i);
						for(idx_t j = vals.Begin(); j != vals.End(); j = vals.Next(j) ) if( j != i )
							rowsum -= vals.Get(j);
						rownorm = 0;
						for(idx_t j = vals.Begin(); j != vals.End(); j = vals.Next(j) )
							rownorm += std::fabs(vals.Get(j));
						if( std::fabs(rowsum) > 1.0e-6*rownorm ) 
							std::cout << __FILE__ << ":" << __LINE__ << " Error: nonequal row sum " << i << " sum " << rowsum << std::endl;
					}
					if( vals.Get(i) )
					{
						
						for(idx_t j = vals.Begin(); j != vals.End(); j = vals.Next(j) ) if( j != i )
						{
							wgt = -vals.Get(j)/vals.Get(i);
							if( wgt ) I.PushBack(id[j],wgt);
						}
					}
					else 
						std::cout << __FILE__ << ":" << __LINE__ << " Error: no diagonal at " << i << std::endl;
					
					I.FinalizeRow();
					
					if( I.RowSize(i) == 0 )
						std::cout << __FILE__ << ":" << __LINE__ << " Error: empty row " << i << std::endl;
					list.Clear();
					vals.Clear();
				}
				else if( CF[i] == C )
				{
					I.PushBack(id[i],1.0);
					I.FinalizeRow();
				}
			}
		}
	}
	void ComputeSizes(const std::vector<char> & CF, idx_t & C_Size, idx_t & F_Size) const
	{
		idx_t S[2] = {0,0};
		for(idx_t i = 0; i < ptr_A->Size(); ++i)
			S[CF[i]-1]++;
		C_Size = S[C-1];
		F_Size = S[F-1];
	}
	bool CheckCols() const
	{
		std::vector<bool> cols(B.Size(),false);
		for(idx_t k = 0; k < B.Size(); ++k)
		{
			for(idx_t j = 0; j < B.RowSize(k); ++j)
				cols[B.Col(k,j)] = true;
		}
		for(idx_t k = 0; k < B.Size(); ++k)
			if( !cols[k] ) 
			{
				std::cout << "no column " << k << std::endl;
				return false;
			}
		return true;
	}
	bool Setup(const CSRMatrix & A)
	{
		bool   print        = GetParameters().Get<int>("verbosity") ? true : false;
		bool   check        = GetParameters().Get<int>("check") ? true : false;
		bool   order        = GetParameters().Get<int>("order") ? true : false;
		bool   write_matrix = GetParameters().Get<int>("write_matrix") ? true : false;
		double phi          = GetParameters().Get<double>("phi");
		int    level        = GetParameters().Get<int>("level");
		idx_t FSize;
		ptr_A = &A;
		if( print ) std::cout << "level " << level << std::endl;
		//compute CF-splitting
		{
			CSRMatrixType<idx_t> S;
			std::vector<char> CF;
			//extract strong graph
			{
				std::vector<double> diag(A.Size(),0.0);
				A.Diagonal(diag);
				StrongGraph(A,diag,phi,S);
			}
			CFSplitting(S,S.Transpose(true),CF);
			ComputeSizes(CF,CSize,FSize);
			if( print ) std::cout << "C size " << CSize << " F size " << FSize << " A size " << A.Size() << "." << std::endl;
			CFRefine(A,S,CF);
			ComputeSizes(CF,CSize,FSize);
			if( print ) std::cout << "C size " << CSize << " F size " << FSize << " A size " << A.Size() << " (refine)." << std::endl;
			if( print ) std::cout << "Compute operators." << std::endl;
			if( print ) std::cout << "Prolongation." << std::endl;
			Operator(A,S,CF,P);
		}
		if( print ) std::cout << "Restriction." << std::endl;
		R = P.Transpose();
		
		G.GetParameters() = GetParameters().SubParameters("Smoother");
		G.Setup(A);
		
		
		if( CSize > 4*FSize || CSize < 16 )// || level > 9 ) 
		{
			if( print ) std::cout << "Setup last level." << std::endl;
			Next = new CoarsestSolver();
			Next->GetParameters() = GetParameters().SubParameters("CoarsestSolver");
			Next->GetParameters().Set("last",1);
		}
		else
		{
			if( print ) std::cout << "Enter level " << level << "." << std::endl;
			Next = new AMGRugeStuben();
			Next->GetParameters() = GetParameters();
			Next->GetParameters().SetRecursive("level",level+1);
			Next->GetParameters().Set("last",0);
		}
		if( print ) std::cout << "Compute coarse system." << std::endl;
		B = R*A*P;
		if( print ) std::cout << "Coarse system size, rows: " << B.Size() << " cols: " << B.Columns() << "." << std::endl;
		if( order )
		{
			if( print ) std::cout << "Sort rows in the coarse system." << std::endl;
			B.SortRows();
		}
		if( check ) 
			CheckCols();
		if( write_matrix )
		{
			if( print ) std::cout << "save A" << level+1 <<".mtx" << std::endl;
			B.Save("A"+to_string(level+1)+".mtx");
			if( print ) std::cout << "save P" << level <<".mtx" << std::endl;
			P.Save("P"+to_string(level)+".mtx");
			if( print ) std::cout << "save R" << level <<".mtx" << std::endl;
			R.Save("R"+to_string(level)+".mtx");
		}
		if( print ) std::cout << "Setup coarse system." << std::endl;
		bool success = Next->Setup(B); 
		rn.resize(CSize);
		xn.resize(CSize);
		r.resize(A.Size());
		return success;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		std::string cycle = GetParameters().Get<std::string>("cycle");
		const CSRMatrix & A = *ptr_A;
		bool success = true;
		std::fill(xn.begin(),xn.end(),0.0);
		//V-cycle
		G.Solve(b,x); // (pre-smoothing)
		std::copy(b.begin(),b.end(),r.begin());
		A.Multiply(-1.0,x,1.0,r); //~ r = b - A*x;
		R.Multiply(1.0,r,0.0,rn); //~ rn = R*(b - A*x);
		success &= Next->Solve(rn,xn);
		P.Multiply(1.0,xn,1.0,x); //~ x= x + P*xn;
		G.Solve(b,x); // (post-(re-)smoothing)
		//V-cycle end
		//W-cycle
		if( cycle == "W" || cycle == "F" )
		{	
			std::copy(b.begin(),b.end(),r.begin());
			A.Multiply(-1.0,x,1.0,r); //~ r = b - A*x;
			R.Multiply(1.0,r,0.0,rn); //~ rn = R*(b - A*x);
			Zero(xn);
			if( cycle == "F" && Next->GetParameters().Get<int>("last") == 0 ) 
			{
				Next->GetParameters().Set("cycle","V");
				success &= Next->Solve(rn,xn);
				Next->GetParameters().Set("cycle","F");
			}
			else success &= Next->Solve(rn,xn);
			P.Multiply(1.0,xn,1.0,x); //~ x+= P*xn;
			G.Solve(b,x); // (post-smoothing)
		}
		//W-cycle end
		return success;
	}
	size_t Bytes() const {return sizeof(const CSRMatrix *) + G.Bytes() + B.Bytes() + P.Bytes() + R.Bytes() + Next->Bytes() + get_bytes(rn) + get_bytes(xn) + get_bytes(r);}
};

#endif //_AMG_RUGE_STUBEN_H
