#ifndef _WRCM_H
#define _WRCM_H
#include "method.h"
#include "priority_queue.h"
#include "graph.h"
#include <cmath>
#include <deque>
/*
 * WeightedReverseCuthillMckee
 * implemented following
 */
struct WRCM_Comparator
{
	std::vector<double>& weights;
public:
	WRCM_Comparator(std::vector<double>& weights) : weights(weights) {}
	WRCM_Comparator(const WRCM_Comparator& b) : weights(b.weights) {}
	WRCM_Comparator& operator = (WRCM_Comparator const& b) { weights = b.weights; return *this; }
	bool operator () (idx_t i, idx_t j) {return weights[i] < weights[j];}
};

template<typename Solver>
class WeightedReverseCuthillMckee : public Methods
{
	std::vector<idx_t> P;
	mutable std::vector<double> Pb, y;
	CSRMatrix B;
	Solver S;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","WeightedReverseCuthillMckee");
		ret.Set("level","*");
		ret.Set("verbosity",1);
		ret.Set("write_matrix",0);
		ret.SubParameters("Solver") = Solver::DefaultParameters();
		return ret;
	}
	WeightedReverseCuthillMckee() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A)
	{
		bool print        = GetParameters().template Get<int>("verbosity")? true : false;
		bool write_matrix = GetParameters().template Get<int>("write_matrix") ? true : false;
		int level         = GetParameters().template Get<int>("level");
		{
			CSRGraph G(A.get_ia(), A.get_ja());
			std::vector<double> weights(A.Size());

			G = G + G.Transpose();

			for (idx_t k = 0; k < A.Size(); ++k)
			{
				for (idx_t jt = 0; jt < A.RowSize(k); ++jt)
				{
					idx_t j = A.Col(k, jt);
					if (k != j)
					{
						weights[k] += fabs(A.Val(k, jt));
						weights[j] += fabs(A.Val(k, jt));
					}
				}
			}
			idx_t index = 0;
			idx_t cur;
			std::deque<idx_t> q;
			std::vector<idx_t> conns;
			P.resize(A.Size());
			std::fill(P.begin(), P.end(), A.Size());
			do
			{
				cur = A.Size();
				for (idx_t k = 0; k < A.Size(); ++k)
				{
					if (P[k] == A.Size())
					{
						cur = k;
						break;
					}
				}
				if (cur == A.Size()) 
					break;
				for (idx_t k = cur + 1; k < A.Size(); ++k)
				{
					if (P[k] == A.Size() && WRCM_Comparator(weights)(k, cur))
						cur = k;
				}
				P[cur] = index++;
				q.push_back(cur);
				while (!q.empty())
				{
					cur = q.front();
					q.pop_front();
					for (idx_t it = 0; it < G.RowSize(cur); ++it)
						if (P[G.Col(cur,it)] == A.Size())
							conns.push_back(G.Col(cur, it));
					std::sort(conns.begin(), conns.end(), WRCM_Comparator(weights));
					for (size_t k = 0; k < conns.size(); ++k)
					{
						P[conns[k]] = index++;
						q.push_back(conns[k]);
					}
					conns.clear();
				}

			}  while (index < A.Size());

			//reverse
			for (idx_t k = 0; k < A.Size(); ++k)
				P[k] = A.Size() - P[k] - 1;
		}

		//assemble new matrix
		{
			std::vector<idx_t> iP(A.Size());
			for (idx_t k = 0; k < A.Size(); ++k) iP[P[k]] = k;
			for (idx_t k = 0; k < A.Size(); ++k)
			{
				for (idx_t jt = 0; jt < A.RowSize(iP[k]); ++jt)
					B.PushBack(P[A.Col(iP[k], jt)], A.Val(iP[k], jt));
				B.FinalizeRow();
			}
		}


		//record result
		if( write_matrix )
		{
			if( level > 0 )
			{
				if( print ) std::cout << "save A_wrcm_"<<level<<".mtx" << std::endl;
				B.Save("A_wrcm_"+to_string(level)+".mtx");
				if( print ) std::cout << "save P"<<level<<".mtx" << std::endl;
				SaveVector("P"+to_string(level)+".mtx",P);
			}
			else 
			{
				if( print ) std::cout << "save A_wrcm.mtx" << std::endl;
				B.Save("A_wrcm.mtx");
				if( print ) std::cout << "save P.mtx" << std::endl;
				SaveVector("P.mtx",P);
			}
		}
		S.SetParameters(GetParameters().SubParameters("Solver"));
		Pb.resize(A.Size());
		y.resize(A.Size());
		return S.Setup(B);
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
			Pb[P[k]] = b[k];
			y[P[k]] = x[k];
		}
		ret = S.Solve(Pb,y);
		//Q[k] is the new position of unknown k
		for(idx_t k = 0; k < size; ++k)
			x[k] = y[P[k]];
		return ret;
	}
	size_t Bytes() const {return S.Bytes() + B.Bytes() + get_bytes(P);}
};

#endif //_MAXIMAL_TRANSVERSAL_H
