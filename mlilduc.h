#ifndef _ILDUC_H
#define _ILDUC_H


#include "csrmatrix.h"
#include "csrtriangular.h"
#include "csctraversal.h"
#include "row_accumulator.h"
#include "invnorm.h"
#include "method.h"

/*
 * MLILDUC - multi-level inverse-based Crout incomplete LDU factorization
 * Implemented following 
 * [1] "Crout Versions of ILU for General Sparse Matrices"
 * by N.Li, Y.Saad, E.Show
 * [2] "A robust and efficient ILU that incorporates the growth of the inverse triangular factors"
 * by M.Bollhoefer
 * [3] "Multilevel preconditioners constructed from inverse-based ILUs"
 * by M.Bollhoefer, Y. Saad
 */
 
template< template<class> class Preprocessor >
class MLILDUC : public Methods
{
	CSRTriangular L, U; ///< factors
	CSRMatrix E,F,S; ///< parts of the matrix and schur complement
	std::vector<double> D; ///< diagonal
	std::vector<idx_t> iP; ///< reordering
	mutable std::vector<double> f, g, y; ///< work vectors
	idx_t BSize, CSize;
	Methods * Next; ///< next-level solver
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","MLILDUC");
		ret.Set("drop_tolerance",0.01);
		ret.Set("diagonal_tolerance",1.0e-7);
		ret.Set("diagonal_perturbation",1.0e-9);
		ret.Set("pivot_condition",2.5);
		ret.Set("write_matrix",0);
		ret.Set("verbosity",1);
		ret.Set("inverse_estimation",1);
		ret.Set("premature_dropping",1);
		ret.Set("check",0);
		ret.Set("level","*");
		ret.SubParameters("Preprocessor") = Preprocessor<DummySolver>::DefaultParameters();
		return ret;
	}
	MLILDUC() : L(CSRTriangular::LowerCSC), U(CSRTriangular::UpperCSR), Next(NULL) {GetParameters() = DefaultParameters();}
	~MLILDUC() {if(Next) delete Next;}
	bool Setup(const CSRMatrix & Ain)
	{
		bool   print = GetParameters().Get<int>("verbosity") & 1 ? true : false;
		//bool   check = GetParameters().Get<int>("check") ? true : false;
		bool   invest = GetParameters().Get<int>("inverse_estimation") ? true : false;
		bool   predrop = GetParameters().Get<int>("premature_dropping") ? true : false;
		bool   write_matrix = GetParameters().Get<int>("write_matrix") ? true : false;
		double tau = GetParameters().Get<double>("drop_tolerance");
		double kappa = GetParameters().Get<double>("pivot_condition");
		double pert = GetParameters().Get<double>("diagonal_perturbation");
		double dtol = GetParameters().Get<double>("diagonal_tolerance");
		int    level        = GetParameters().Get<int>("level");
		std::vector<bool> pivot(Ain.Size(),false);
		idx_t swaps = 0;
		{
			idx_t  report_pace = std::max<idx_t>(1,Ain.Size() / 25);
			CSRMatrix A = Ain; //have to copy and sort rows due to CSCTraversal
			A.SortRows();
			L.Clear();
			U.Clear();
			CSCTraversal At(A,A.Size()), Lt(L,A.Size()), Ut(U,A.Size());
			RowAccumulator<double> u(A.Size()), l(A.Size());
			std::vector<double> unorms(A.Size(),0.0), lnorms(A.Size(),0.0);
			Invnorm iLest(invest || print ? A.Size() : 0), iUest(invest || print ? A.Size() : 0);
			double iLnorm = 1, iUnorm = 1, Dmax = 0, Dmin = 1.0e+20;
			D.resize(A.Size(),0.0);
			if( print ) 
				std::cout << "tau: " << tau << " inverse estimation: " << (invest ? "yes":"no") << std::endl;
			for (idx_t k = 0; k < A.Size(); k++)
			{
				if( print )
				{
					if( k % report_pace == 0 )
					{
						std::cout << "precond: " << std::setw(12) << (k+1)*100.0 / A.Size();
						std::cout << " L " << std::setw(12) << iLnorm;
						std::cout << " D " << std::setw(12) << Dmax/Dmin;
						std::cout << " U " << std::setw(12) << iUnorm;
						std::cout << " swaps " << std::setw(6) << swaps;
						std::cout << "\r";
						std::cout.flush();
					}
				}
				{
					//Compute U-part
					{
						//Uncompress row
						idx_t curr = k;
						u[k] = 0.0;
						for(idx_t j = 0; j < A.RowSize(k); ++j) if( A.Col(k,j) >= k )
							curr = u.InsertOrdered(curr,A.Col(k,j),A.Val(k,j));
						if( pert )
							u.Get(k) = u.Get(k) * (1.0 + pert) + (u.Get(k) < 0.0 ? -1.0 : 1.0) * pert;
						//U part elimination with L
						for(idx_t j = Lt.Begin(); j != Lt.End(); j = Lt.Next(j))
						{
							double v = -L.Val(j,Lt.Position(j)) * D[j];
							if( !predrop || std::fabs(v)*unorms[j] > tau )
							{
								curr = k;
								for(idx_t r = 0; r < U.RowSize(j); ++r) if( U.Col(j,r) >= k )
									curr = u.InsertOrdered(curr,U.Col(j,r),v*U.Val(j,r));
							}
						}
						if( invest || print )
							iUnorm = iUest.Estimate(u,k);
					}
					//Compute L-part
					{
						//Uncompress column
						idx_t curr = k;
						l[k] = 0.0;
						for(idx_t j = At.Begin(); j != At.End(); j = At.Next(j)) if( j >= k )
							curr = l.InsertOrdered(curr,j,A.Val(j,At.Position(j)));
						At.NextColumn();
						if( pert )
							l.Get(k) = l.Get(k) * (1.0 + pert) + (l.Get(k) < 0.0 ? -1.0 : 1.0) * pert;
						//L part elimination with U
						for(idx_t j = Ut.Begin(); j != Ut.End(); j = Ut.Next(j))
						{
							double v = -U.Val(j,Ut.Position(j)) * D[j];
							if( !predrop || std::fabs(v)*lnorms[j] > tau )
							{
								curr = k;
								for(idx_t r = 0; r < L.RowSize(j); ++r) if( L.Col(j,r) >= k )
									curr = l.InsertOrdered(curr,L.Col(j,r),v*L.Val(j,r));
							}
						}
						if( invest || print )
							iLnorm = iLest.Estimate(l,k);
					}
				}
				if( iUnorm < kappa && iLnorm < kappa )
				{
					//retrive diagonal
					D[k] = (u.Get(k) + l.Get(k))*0.5;
					D[k] = (D[k] < 0.0 ? -1.0 : 1.0)*std::max(std::fabs(D[k]),dtol);
					Dmax = std::max(Dmax,std::fabs(D[k]));
					Dmin = std::min(Dmin,std::fabs(D[k]));
					{
						//Assemble U-part
						{
							//Scale u by diagonal
							double unorm = 0;
							for(idx_t j = u.Begin(); j != u.End(); j = u.Next(j))
							{
								u.Get(j) /= D[k];
								unorm += std::pow(u.Get(j),2);
							}
							unorm = std::sqrt(unorm);
							unorms[k] = unorm;
							//inverse estimation
							if( invest || print )
							{
								iUest.Update(u,k);
								if( invest ) unorm /= iUnorm;
							}
							//assemble row of U
							U.PushBack(k,u.Get(k));
							for(idx_t j = u.Next(u.Begin()); j != u.End(); j = u.Next(j))
							{
								double v = u.Get(j);
								if( std::fabs(v) > tau*unorm ) U.PushBack(j,v);
							}
						}
						//Assemble L-part
						{
							//Scale u by diagonal
							double lnorm = 0;
							for(idx_t j = l.Begin(); j != l.End(); j = l.Next(j))
							{
								l.Get(j) /= D[k];
								lnorm += std::pow(l.Get(j),2);
							}
							lnorm = std::sqrt(lnorm);
							lnorms[k] = lnorm;
							//inverse estimation
							if( invest || print )
							{
								iLest.Update(l,k);
								if( invest ) lnorm /= iLnorm;
							}
							//assemble row of U
							L.PushBack(k,l.Get(k));
							for(idx_t j = l.Next(l.Begin()); j != l.End(); j = l.Next(j))
							{
								double v = l.Get(j);
								if( std::fabs(v) > tau*lnorm ) L.PushBack(j,v);
							}
						}
					}
				}
				else
				{
					pivot[k] = true;
					swaps++;
				}
				{ //Finalize U-part
					U.FinalizeRow();
					Ut.NewRow(k);
					Ut.NextColumn();
					u.Clear();
				}
				{ //Finalize L-part
					L.FinalizeRow();
					Lt.NewRow(k);
					Lt.NextColumn();
					l.Clear();
				}
			}
			if( print )
			{
				std::cout << "      nonzeros in A: " << A.Nonzeros() << " consumed: " << A.Bytes()/1024.0/1024.0 << "Mb " << std::endl;
				std::cout << "      nonzeros in L: " << L.Nonzeros() << " consumed: " << L.Bytes()/1024.0/1024.0 << "Mb " << std::endl;
				std::cout << "      nonzeros in U: " << U.Nonzeros() << " consumed: " << U.Bytes()/1024.0/1024.0 << "Mb " << std::endl;
				std::cout << "      consumed by D: " << D.capacity()*sizeof(double)/1024.0/1024.0 << "Mb " << std::endl;
				size_t bytes = A.Bytes() + L.Bytes() + U.Bytes();
				bytes += (unorms.capacity()+lnorms.capacity()+D.capacity())*sizeof(double);
				bytes += Lt.Bytes() + Ut.Bytes();
				bytes += iUest.Bytes() + iLest.Bytes();
				std::cout << "      consumed in factorization " << bytes/1024.0/1024.0 << "Mb " << std::endl;
				std::cout << "      fill-in LU: " << (L.Nonzeros() + U.Nonzeros())/(1.0*A.Nonzeros()) << std::endl;
				std::cout << "      estimated inverse norms L " << iLnorm << " D " << Dmax/Dmin << " U " << iUnorm << std::endl;
				std::cout << "      swaps " << swaps << std::endl;
			}
		}
		bool success = true;
		if( swaps )
		{
			const CSRMatrix & A = Ain;
			CSize = swaps, BSize = Ain.Size() - CSize;
			std::vector<idx_t> P(A.Size());
			//reorder system
			{
				idx_t id = 0;
				if( print ) std::cout << "Compute ordering, P and iP" << std::endl;
				for(idx_t k = 0; k < A.Size(); ++k) if(!pivot[k] ) P[k] = id++;
				assert(id == BSize);
				for(idx_t k = 0; k < A.Size(); ++k) if( pivot[k] ) P[k] = id++;
				assert(id == A.Size());
				iP.resize(A.Size());
				for(idx_t k = 0; k < A.Size(); ++k) iP[P[k]] = k;
				
				if (print) std::cout << "reorder L" << std::endl;
				//reorder L
				{
					//remove gaps in ia array
					L.RemoveEmptyRows();
					assert(L.Size() == BSize);
					//change column indices
					for(idx_t k = 0; k < L.Size(); ++k)
					{
						for(idx_t j = 0; j < L.RowSize(k); ++j)
							L.Col(k,j) = P[L.Col(k,j)];
					}
					L.SortRows(); //TODO: maybe only move diagonal
					//filter out entries that are outside range
					L.ChopColumns(0,BSize);
				}
				if (print) std::cout << "reorder U" << std::endl;
				//reorder U
				{
					//remove gaps in ia array
					U.RemoveEmptyRows();
					assert(U.Size() == BSize);
					//change row indices
					for(idx_t k = 0; k < U.Size(); ++k)
					{
						for(idx_t j = 0; j < U.RowSize(k); ++j)
							U.Col(k,j) = P[U.Col(k,j)];
					}
					U.SortRows(); //TODO: maybe only move diagonal
					//filter out entries that are outside range
					U.ChopColumns(0,BSize);
				}
				if( print ) std::cout << "reorder D" << std::endl;
				//reorder D
				{
					//compress D
					idx_t id = 0;
					for(idx_t k = 0; k < A.Size(); ++k) 
						if( !pivot[k] )
							D[id++] = D[k];
					assert(id == BSize);
					D.resize(id);
				}
				if( print ) std::cout << "assemble blocks" << std::endl;
				{
					/*
					 *     | B   F |
					 * A = |       |
					 *     | E   C |
					 * 
					 * S = C - E B^{-1} F
					 */
					if( print ) std::cout << "assemble F" << std::endl;
					{//assemble F (upper-right block) in row-major format
						for(idx_t k = 0; k < BSize; ++k)
						{
							idx_t iPk = iP[k];
							for(idx_t j = 0; j < A.RowSize(iPk); ++j)
							{
								idx_t Pc = P[A.Col(iPk,j)];
								if( Pc >= BSize ) 
									F.PushBack(Pc-BSize,A.Val(iPk,j));
							}
							F.FinalizeRow();
						}
						F.SortRows();
					}
					if( print ) std::cout << "assemble E" << std::endl;
					{ //assemble E (lower-left block) in row-major format
						for(idx_t k = BSize; k < A.Size(); ++k)
						{
							idx_t iPk = iP[k];
							for(idx_t j = 0; j < A.RowSize(iPk); ++j)
							{
								idx_t Pc = P[A.Col(iPk,j)];
								if( Pc < BSize ) 
									E.PushBack(Pc,A.Val(iPk,j));
							}
							E.FinalizeRow();
						}
						E.SortRows();
					}
				}
			}
			//compute Schur, S-version
			{
				CSRMatrix EU, LF;
				{ //compute EU
					if( print ) std::cout << "solve for EU^{-1}" << std::endl;
					{
						RowAccumulator<double> row(BSize);
						for(idx_t k = 0; k < CSize; ++k)
						{
							//unpack row of E
							idx_t curr = 0;
							for(idx_t j = 0; j < E.RowSize(k); ++j)
								curr = row.InsertOrdered(curr,E.Col(k,j),E.Val(k,j)/U.Val(E.Col(k,j),0));
							//elimination with U
							for(idx_t j = row.Begin(); j != row.End(); j = row.Next(j))
							{
								if( 1 + row.Get(j) != 1 )
								{
									curr = j;
									assert(U.RowSize(j) != 0); //at least diagonal
									assert(U.Col(j,0) == j); //diagonal goes first
									for(idx_t l = 1; l < U.RowSize(j); ++l)
									{
										assert(U.Col(U.Col(j,l),0) == U.Col(j,l)); //diagonal goes first
										curr = row.InsertOrdered(curr,U.Col(j,l),-row.Get(j)*U.Val(j,l)/U.Val(U.Col(j,l),0));
									}
								}
							}
							row.Get(EU.get_ja(),EU.get_a());
							row.Clear();
							EU.FinalizeRow();
						}
						assert(EU.Size() == CSize);
					}
					if( print ) std::cout << "Size: " << EU.Size() << " Nonzeroes: " << EU.Nonzeros() << std::endl;
				}
				{ //compute LF
					if( print ) std::cout << "solve for L^{-1}F" << std::endl;
					//eliminate rows of LF by L column
					{
						CSCTraversal Lt(L, BSize);
						RowAccumulator<double> row(CSize);
						for (idx_t k = 0; k < BSize; ++k)
						{
							//fill LF row with F row
							idx_t curr = 0;
							for (idx_t j = 0; j < F.RowSize(k); ++j)
								curr = row.InsertOrdered(curr, F.Col(k, j), F.Val(k, j));
							//eliminate LF rows from current row
							for (idx_t j = Lt.Begin(); j != Lt.End(); j = Lt.Next(j)) //until diagonal
							{
								if (j != k)
								{
									curr = 0;
									for (idx_t l = 0; l < LF.RowSize(j); ++l)
										curr = row.InsertOrdered(curr, LF.Col(j, l), -LF.Val(j, l) * L.Val(j, Lt.Position(j)));
								}
							}
							for (idx_t j = row.Begin(); j != row.End(); j = row.Next(j))
								row.Get(j) /= L.Val(k, 0);
							row.Get(LF.get_ja(), LF.get_a());
							row.Clear();
							LF.FinalizeRow();
							Lt.NextColumn();
						}
					}
					if( print ) std::cout << "Size: " << LF.Size() << " Nonzeroes: " << LF.Nonzeros() << std::endl;
				}
				{//S = C - EU*LF;
					if( print ) std::cout << "Compute Schur" << std::endl;
					RowAccumulator<double> row(CSize);
					for(idx_t k = BSize; k < A.Size(); ++k)
					{
						//first assemble C (lower-right block) into row accumulator
						idx_t iPk = iP[k];
						for(idx_t j = 0; j < A.RowSize(iPk); ++j)
						{
							idx_t Pc = P[A.Col(iPk,j)];
							if( Pc >= BSize )
								row[Pc-BSize] = A.Val(iPk,j);
						}
						row.Sort();
						//for k-th row of EU add all rows of LF (scaled by D) to row-accumulator
						for(idx_t jt = 0; jt < EU.RowSize(k-BSize); ++jt)
						{
							idx_t j = EU.Col(k-BSize,jt);
							idx_t curr = 0;
							for (idx_t lt = 0; lt < LF.RowSize(j); ++lt)
								curr = row.InsertOrdered(curr, LF.Col(j, lt), -EU.Val(k - BSize, jt) * LF.Val(j, lt) / D[j]);
						}
						//put row accumulator to S
						row.Get(S.get_ja(),S.get_a());
						row.Clear();
						S.FinalizeRow();
					}
					assert(S.Size() == CSize);
					if( print ) std::cout << "Size: " << S.Size() << " Nonzeroes: " << S.Nonzeros() << std::endl;
				}
			}
			//setup next level system
			if( print ) std::cout << "Setup next level" << std::endl;
			Next = new Preprocessor<MLILDUC>();
			Next->GetParameters() = GetParameters();
			Next->GetParameters() = GetParameters().SubParameters("Preprocessor");
			Next->GetParameters().SubParametersSearchRecursive("Solver", "name", "DummySolver") = GetParameters();
			Next->GetParameters().SetRecursive("level",level+1);
			success = Next->Setup(S);
			//allocate vectors
			//if( print ) std::cout << "Allocate vectors" << std::endl;
			f.resize(BSize);
			g.resize(CSize);
			y.resize(CSize);
		}
		if( write_matrix )
		{
			if( print ) std::cout << "save L" << level << ".mtx" << std::endl;
			L.Save("L"+to_string(level)+".mtx");
			if( print ) std::cout << "save U" << level << ".mtx" << std::endl;
			U.Save("U"+to_string(level)+".mtx");
			if( print ) std::cout << "save E" << level << ".mtx" << std::endl;
			E.Save("E"+to_string(level)+".mtx");
			if( print ) std::cout << "save F" << level << ".mtx" << std::endl;
			F.Save("F"+to_string(level)+".mtx");
			if( print ) std::cout << "save S" << level << ".mtx" << std::endl;
			S.Save("L"+to_string(level)+".mtx");
			if (print) std::cout << "save D" << level << ".mtx" << std::endl;
			SaveVector("D" + to_string(level) + ".mtx",D);
			if (print) std::cout << "save iP" << level << ".mtx" << std::endl;
			SaveVector("iP" + to_string(level) + ".mtx", iP);
		}
		return  success;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		//bool print = GetParameters().Get<int>("verbosity") & 2 ? true : false;
		if( Next )
		{
			//the original matrix A was separated by multilevel algorithm to the following form
			//     | B  F |
			// A = |      |
			//     | E  C |
			// In order to apply solve phase on this matrix, we consider the matrix in following sense:
			//     | I         0 | | B   0 | | I    B^{-1} F |
			// A = |             | |       | |               | = L D U
			//     | E B^{-1}  I | | 0   S | | 0           I |
			// where S = C - E B^{-1} F
			// consider the system
			//  | B  F | | u |   | f |
			//  |      | |   | = |   |
			//  | E  C | | y |   | g |
			// then solution is obtained by steps:
			// 1) ~f = B^{-1} f
			// 2) ~g = g - E ~f
			// 3) y = S^{-1} ~g
			// 4) u = ~f - B^{-1} F y
			for(idx_t k = 0; k < BSize; ++k) f[k] = b[iP[k]];
			for(idx_t k = BSize; k < BSize+CSize; ++k) g[k-BSize] = b[iP[k]];
			L.Solve(f);
			for(idx_t k = 0; k < BSize; ++k) f[k] /= D[k];
			U.Solve(f);
			for(idx_t k = 0; k < BSize; ++k) x[iP[k]] = f[k];
			E.Multiply(-1,f,1,g);
			Next->Solve(g,y);
			for(idx_t k = BSize; k < BSize+CSize; ++k) x[iP[k]] = y[k-BSize];
			F.Multiply(1,y,0,f);
			L.Solve(f);
			for(idx_t k = 0; k < BSize; ++k) f[k] /= D[k];
			U.Solve(f);
			for(idx_t k = 0; k < BSize; ++k) x[iP[k]] -= f[k];
		}
		else
		{
			x.resize(b.size());
			std::copy(b.begin(),b.end(),x.begin());
			L.Solve(x);
			for(idx_t k = 0; k < D.size(); ++k) x[k] /= D[k];
			U.Solve(x);
		}
		return true;
	}
	size_t Bytes() const {return U.Bytes() + L.Bytes() + get_bytes(D) + sizeof(const CSRMatrix *);}
};

#endif //_ILDUC_H