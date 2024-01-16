#ifndef _MAXIMAL_TRANSVERSAL_H
#define _MAXIMAL_TRANSVERSAL_H
#include "method.h"
#include "priority_queue.h"
#include <cmath>
/*
 * MaximalTransversal
 * implemented following
 * "A New Pivoting Strategy for Gaussian Elimination"
 * by M.Olschowka and A.Neumaier
 * and
 * "Design, Implementation and Analysis of Maximum Transversal Algorithms"
 * by I.S.Duff, K.Kaya and B.Ucar
 */
template<typename Solver>
class MaximalTransversal : public Methods
{
	std::vector<idx_t> Q;
	std::vector<double> L, R;
	mutable std::vector<double> Lb, y;
	CSRMatrix B;
	Solver S;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","MaximalTransversal");
		ret.Set("maxiters",5);
		ret.Set("level","*");
		ret.Set("verbosity",1);
		ret.Set("write_matrix",0);
		ret.Set("reorder",1);
		ret.Set("check",1);
		ret.SubParameters("Solver") = Solver::DefaultParameters();
		return ret;
	}
	MaximalTransversal() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		bool print        = GetParameters().template Get<int>("verbosity")? true : false;
		bool write_matrix = GetParameters().template Get<int>("write_matrix") ? true : false;
		bool reorder      = GetParameters().template Get<int>("reorder")? true : false;
		bool check        = GetParameters().template Get<int>("check")? true : false;
		int maxiters      = GetParameters().template Get<int>("maxiters");
		int level         = GetParameters().template Get<int>("level");
		const idx_t NDEF = std::numeric_limits<idx_t>::max(), EOL = NDEF-1;
		idx_t size = A.Size();
		PriorityQueue<double, std::less_equal<double> > queue(size,std::numeric_limits<double>::max());
		std::vector<idx_t> perm(size,NDEF), iperm(size,NDEF), parent(size,NDEF), clist(size,NDEF), gaps;
		std::vector< std::pair<idx_t,idx_t> > col_pos(size), aug_pos(size);
		std::vector<double> Bmax(size,0.0);
		std::vector<double> U(size,std::numeric_limits<double>::max());
		std::vector<double> V(size,std::numeric_limits<double>::max());
		double min_path, aug_path;
		idx_t end_path, trace, clist_beg;
		idx_t k, i, j, m, l, p;
		B = A;
		// Initial LOG transformation to dual problem and initial extreme match
		for(idx_t k = 0; k < size; ++k)
		{
			for(idx_t j = 0; j < B.RowSize(k); ++j)
			{
				idx_t i = B.Col(k,j);
				B.Val(k,j) = fabs(B.Val(k,j));
				Bmax[i] = std::max(Bmax[i],B.Val(k,j));
			}
		}
		for(idx_t k = 0; k < size; ++k)
		{
			for (idx_t j = 0; j < B.RowSize(k); ++j)
			{
				idx_t i = B.Col(k,j);
				if( Bmax[i] == 0 || B.Val(k,j) == 0 )
					B.Val(k,j) = std::numeric_limits<double>::max();
				else
				{
					B.Val(k,j) = log(Bmax[i])-log(B.Val(k,j));
					U[i] = std::min(U[i],B.Val(k,j));
				}
			}
		}
		for(idx_t k = 0; k < size; ++k)
		{
			for (idx_t j = 0; j < B.RowSize(k); ++j)
			{
				idx_t i = B.Col(k,j);
				double u = B.Val(k,j) - U[i];
				V[k] = std::min(V[k],u);
			}
		}
		// Update cost and match
		for(idx_t k = 0; k < size; ++k)
		{
			for (idx_t j = 0; j < B.RowSize(k); ++j)
			{
				idx_t i = B.Col(k,j);
				double u = fabs(B.Val(k,j) - V[k] - U[i]);
				if( u < 1.0e-30 )
				{
					if( perm[i] == NDEF )
					{
						perm[i] = k;
						iperm[k] = i;
						col_pos[k] = std::make_pair(k,j);
						break;
					}
				}
			}
		}
		// 1-step augmentation
		for(idx_t k = 0; k < size; ++k)
		{
			if( iperm[k] == NDEF ) //unmatched row
			{
				for (idx_t j = 0; j < B.RowSize(k); ++j)
				{
					idx_t i = B.Col(k,j);
					double u = fabs(B.Val(k,j) - V[k] - U[i]);
					if( u <= 1.0e-30 )
					{
						idx_t m = perm[i];
						if( m == NDEF ) throw "Unexpected value";
						// Search other row in C for 0
						for(idx_t l = 0; l < B.RowSize(m); ++l)
						{
							idx_t p = B.Col(m,l);
							double u = fabs(B.Val(m,l) - V[m] - U[p]);
							if( u <= 1.0e-30 )
							{
								if( perm[p] == NDEF )
								{
									perm[i] = k;
									iperm[k] = i;
									col_pos[k] = std::make_pair(k,j);
									perm[p] = m;
									iperm[m] = p;
									col_pos[m] = std::make_pair(m,l);
									break;
								}
							}
						}
						if( iperm[k] != NDEF ) break;
					}
				}
			}
		}
		// Weighted bipartite matching
		for(k = 0; k < size; ++k)
		{
			if( iperm[k] != NDEF ) 
				continue;
			l = k;
			clist_beg = EOL;
			parent[l] = NDEF;
			end_path = NDEF;
			trace = k;
			min_path = 0;
			aug_path = std::numeric_limits<double>::max();
			// fill queue
			while(true)
			{
				for(j = 0; j < B.RowSize(l); ++j)
				{
					i = B.Col(l,j);
					if( clist[i] != NDEF ) continue;
					double len = fabs(min_path + B.Val(l,j) - V[l] - U[i]);
					if( len < aug_path )
					{
						if( perm[i] == NDEF )
						{
							end_path = i;
							trace = l;
							aug_path = len;
							aug_pos[i] = std::make_pair(l,j);
						}
						else if( len < queue.GetKey(i) )
						{
							parent[perm[i]] = l;
							aug_pos[i] = std::make_pair(l,j);
							if( queue.Contains(i) )
								queue.DecreaseKey(i,len);
							else
								queue.Push(i,len);
						}
					}
				}
				if( queue.Empty() ) break;
				p = queue.Pop();
				min_path = queue.GetKey(p);
				if( aug_path <= min_path ) 
				{
					queue.ResetKey(p,std::numeric_limits<double>::max());
					break;
				}
				clist[p] = clist_beg;
				clist_beg = p;
				l = perm[p];
			}
			if( end_path != NDEF )
			{
				i = clist_beg;
				while(i != EOL)
				{
					U[i] += queue.GetKey(i) - aug_path;
					if( perm[i] != NDEF ) 
						V[perm[i]] = B.Val(col_pos[perm[i]]) - U[i];
					queue.ResetKey(i,std::numeric_limits<double>::max());
					//clear list
					m = clist[i];
					clist[i] = NDEF;
					i = m;
				}

				i = end_path;
				while(trace != NDEF)
				{
					m = iperm[trace];
					perm[i] = trace;
					iperm[trace] = i;
					col_pos[trace] = aug_pos[i];
					V[trace] = B.Val(col_pos[trace]) - U[i];
					i = m;
					trace = parent[trace];
				}
				while( !queue.Empty() )
					queue.ResetKey(queue.Pop(),std::numeric_limits<double>::max());
			}
		}		
		//fix gaps in permutation
		for(idx_t k = 0; k < size; ++k) 
			iperm[k] = NDEF;
		for(idx_t k = 0; k < size; ++k)
		{
			if( perm[k] != NDEF )
				iperm[perm[k]] = 0;
		}
		for(idx_t k = 0; k < size; ++k)
			if( iperm[k] == NDEF )
				gaps.push_back(k);
		if( print )
			std::cout << "gaps: " << gaps.size() << std::endl;
		for(idx_t k = 0; k < size; ++k)
		{
			if( perm[k] == NDEF )
			{
				perm[k] = gaps.back();
				gaps.pop_back();
			}
		}
		//check correctness
		if( check )
		{
			for (idx_t k = 0; k < size; ++k)
			{
				double v = 1;
				if( V[k] != std::numeric_limits<double>::max() ) 
					v = exp(V[k]);
				//flip sign and check rescaling
				for (idx_t j = 0; j < A.RowSize(k); ++j)
				{
					idx_t i = A.Col(k,j);
					double u = 1;
					if( U[i] != std::numeric_limits<double>::max() && Bmax[i] != 0 ) 
						u = exp(U[i])/Bmax[i];
					if( fabs(v*A.Val(k,j)*u) > 1 + 1.0e-7 )
					{
						std::cout << "element(" << k << "," << i << ") " << A.Val(k,j);
						std::cout << " scaled " << v*A.Val(k,j)*u;
						std::cout << " u " << u << " v " << v;
						std::cout << " U " << U[i] << " V " << V[k] << " Bmax " << Bmax[i];
						std::cout << std::endl;
					}
				}
			}
		}
		//assemble rescaling
		L.resize(size);
		R.resize(size);
		for (idx_t k = 0; k < size; ++k)
		{
			double v = 1, u = 1;
			if( V[k] != std::numeric_limits<double>::max() ) 
				v = exp(V[k]);
			if( U[k] != std::numeric_limits<double>::max() && Bmax[k] != 0 ) 
				u = exp(U[k])/Bmax[k];
			for (idx_t j = 0; j < A.RowSize(k); ++j)
			{
				if( perm[A.Col(k,j)] == k && v*A.Val(k,j)*u < 0.0 ) 
					v *= -1;
			}
			L[k] = v;
			R[k] = u;
		}
		if( !reorder )
		{
			for (k = 0; k < size; ++k)
				perm[k] = k;
		}
		//rescale matrix
		for (idx_t k = 0; k < size; ++k)
		{
			for (idx_t j = 0; j < A.RowSize(k); ++j)
			{
				idx_t i = A.Col(k,j);
				B.Col(k,j) = perm[i];
				B.Val(k,j) = L[k]*A.Val(k,j)*R[i];
			}
		}
		//save permutation
		Q.resize(size);
		idx_t perms = 0;
		for(idx_t k = 0; k < size; ++k)
		{
			Q[k] = perm[k];
			if( Q[k] != k ) 
				perms++;
		}
		if( print ) 
			std::cout << "total permutations " << perms << " / " << size << std::endl;
		//improve I-dominance
		{
			for(idx_t k = 0; k < size; ++k)
			{
				for (idx_t j = 0; j < B.RowSize(k); ++j)
				{
					double u = fabs(B.Val(k,j));
					if( u > 1.0e-12 )
						B.Val(k,j) = -log(u);
					else
						B.Val(k,j) = std::numeric_limits<double>::max();
				}
			}
			for(idx_t k = 0; k < size; ++k) 
				Bmax[k] = 0.0;
			for(int l = 0; l < maxiters; l++)
			{
				for(idx_t k = 0; k < size; ++k) 
					U[k] = V[k] = std::numeric_limits<double>::max();
				for(idx_t k = 0; k < size; ++k)
				{
					for (idx_t j = 0; j < B.RowSize(k); ++j)
					{
						idx_t i = B.Col(k,j); //column number
						if( i != k ) //out of diagonal
						{
							double u = B.Val(k,j) + Bmax[k] - Bmax[i];
							U[k] = std::min(U[k],u);// update Y1
							V[i] = std::min(V[i],u);// update Y2
						}
					}
				}
				for (idx_t k = 0; k < size; ++k)
				{
					if( U[k] != std::numeric_limits<double>::max() &&
						V[k] != std::numeric_limits<double>::max() )
						Bmax[k] += (V[k]-U[k])*0.5;
				}
			}
			for (idx_t k = 0; k < size; ++k)
			{
				U[k] = exp(-Bmax[k]);
				V[k] = exp(Bmax[k]);
				if( U[k] != U[k] ) U[k] = 1;
				if( V[k] != V[k] ) V[k] = 1;
			}
			for (idx_t k = 0; k < size; ++k)
			{
				L[k] *= U[k];
				R[k] *= V[Q[k]];
			}
			for (idx_t k = 0; k < size; ++k)
			{
				for (idx_t j = 0; j < B.RowSize(k); ++j)
					B.Val(k,j) = L[k] * A.Val(k,j) * R[A.Col(k,j)];
			}
		}
		

		//record result
		if( write_matrix )
		{
			if( level > 0 )
			{
				if( print ) std::cout << "save B"<<level<<".mtx" << std::endl;
				B.Save("B"+to_string(level)+".mtx");
				if( print ) std::cout << "save L"<<level<<".mtx" << std::endl;
				SaveVector("L"+to_string(level)+".mtx",L);
				if( print ) std::cout << "save R"<<level<<".mtx" << std::endl;
				SaveVector("R"+to_string(level)+".mtx",R);
				if( print ) std::cout << "save Q"<<level<<".mtx" << std::endl;
				SaveVector("Q"+to_string(level)+".mtx",Q);
			}
			else 
			{
				if( print ) std::cout << "save B.mtx" << std::endl;
				B.Save("B.mtx");
				if( print ) std::cout << "save L.mtx" << std::endl;
				SaveVector("L.mtx",L);
				if( print ) std::cout << "save R.mtx" << std::endl;
				SaveVector("R.mtx",R);
				if( print ) std::cout << "save Q.mtx" << std::endl;
				SaveVector("Q.mtx",Q);
			}
		}
		S.SetParameters(GetParameters().SubParameters("Solver"));
		bool success = S.Setup(B);
		Lb.resize(A.Size());
		y.resize(A.Size());
		return success;
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		idx_t size = B.Size();
		// A x = b
		// L A x = L b
		// (L A R Q) y = L b, x = R Q y
		//
		// R Q y = x
		// Q y = R^{-1} x
		//
		// B y = L b
		// y = B^{-1} L b
		// x = R Q y
		assert(b.size() == size);
		bool ret = false;
		x.resize(size,0.0);
		for(idx_t k = 0; k < size; ++k)
		{
			Lb[k] = L[k]*b[k];
			y[Q[k]] = x[k]/R[k];
		}
		ret = S.Solve(Lb,y);
		//Q[k] is the new position of unknown k
		for(idx_t k = 0; k < size; ++k)
			x[k] = y[Q[k]]*R[k];
		return ret;
	}
	size_t Bytes() const {return S.Bytes() + B.Bytes() + get_bytes(Q) + get_bytes(L) + get_bytes(R) + get_bytes(Lb) + get_bytes(y);}
};

#endif //_MAXIMAL_TRANSVERSAL_H
