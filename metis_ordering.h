#ifndef _METIS_ORDERING_H
#define _METIS_ORDERING_H
#if defined(USE_METIS)
//same as in cmake
#define IDXTYPEWIDTH 32
#define REALTYPEWIDTH 64
#define idx_t  metis_idx_t
#define real_t metis_real_t
#include <metis.h>
#undef idx_t
#undef real_t
#else
#define metis_idx_t  int
#define metis_real_t float
#define METIS_NOPTIONS 40
#endif
#include "method.h"
#include "graph.h"
#include "priority_queue.h"
/*
 * Metis
 * 
 * call to metis partitioner to reorder the system
 * Todo: METIS_NodeNDP
 */
template<typename Solver>
class Metis : public Methods
{
	std::vector<idx_t> Q;
	CSRMatrix B;
	Solver S;
	mutable std::vector<double> Qb, y;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","Metis");
		ret.Set("level","*");
		ret.Set("verbosity",1);
		ret.Set("sorted",1);
		ret.Set("method", "recursive");
		//ret.Set("method", "nodend");
		ret.Set("method_info", "kway,recursive,nodend");
		ret.Set("parts", 4);
		ret.Set("config", 0);
		ret.Set("debug", 0);
		ret.Set("separator", 1);
		ret.Set("write_matrix",1);
		ret.Set("write_format", "mtx");
		ret.SubParameters("Solver") = Solver::DefaultParameters();
		return ret;
	}
	Metis() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix& A)
	{
		bool print = GetParameters().template Get<int>("verbosity") ? true : false;
		bool write_matrix = GetParameters().template Get<int>("write_matrix") ? true : false;
		bool sorted = GetParameters().template Get<int>("sorted") ? true : false;
		int level = GetParameters().template Get<int>("level");
		std::string mt = GetParameters().template Get<std::string>("method");
		metis_idx_t nparts = GetParameters().template Get<metis_idx_t>("parts");
		idx_t size = A.Size();
		//output data - which entity should belong to which part
		std::vector<metis_idx_t> part;
		std::vector<idx_t> blocks;
		{ //metis memory
			//graph assembly
			CSRGraph G(A.get_ia(), A.get_ja());
			if (!G.Symmetric())
				G = G + G.Transpose();
			{ // metis input graph
				//These arrays determine local sparse matrix
				std::vector<metis_idx_t> xadj;
				std::vector<metis_idx_t> adjncy;
				//Number of constraints
				metis_idx_t ncon = 1;
				//set options
				metis_idx_t options[METIS_NOPTIONS];
				//number of edgecuts
				metis_idx_t edgecut = 0;
				//graph size
				metis_idx_t nvtxs = A.Size();
				//for NodeND
				std::vector<metis_idx_t> perm, iperm;

				if (print)
					std::cout << "Assemble graph" << std::endl;

				{ //graph copy
					xadj.resize(G.get_ia().size());
					std::copy(G.get_ia().begin(), G.get_ia().end(), xadj.begin());
					adjncy.resize(G.get_ja().size());
					std::copy(G.get_ja().begin(), G.get_ja().end(), adjncy.begin());
				}

				if (print)
					std::cout << "Call METIS" << std::endl;

#if defined(USE_METIS)
				METIS_SetDefaultOptions(options);
				options[METIS_OPTION_NUMBERING] = 0;
				if (GetParameters().Get<int>("debug")) //from INMOST, maybe wrong for METIS
				{
					options[0] = 1;
					options[1] = 1;
					options[2] = 15;
				}
				options[METIS_OPTION_CONTIG] = GetParameters().Get<metis_idx_t>("config");
				if (nparts == 1)
				{
					part.resize(size, 0);
					for (idx_t k = 0; k < part.size(); ++k)
						part[k] = 0;
				}
				else if (mt == "recursive")
				{
					part.resize(size, 0);
					METIS_PartGraphRecursive(&nvtxs, &ncon, &xadj[0],
						&adjncy[0], NULL, NULL, NULL,
						&nparts, NULL, NULL, options,
						&edgecut, &part[0]);
				}
				else if (mt == "kway")
				{
					part.resize(size, 0);
					METIS_PartGraphKway(&nvtxs, &ncon, &xadj[0],
						&adjncy[0], NULL, NULL, NULL,
						&nparts, NULL, NULL, options,
						&edgecut, &part[0]);
				}
				else if (mt == "nodend")
				{
					perm.resize(size);
					iperm.resize(size);
					METIS_NodeND(&nvtxs, &xadj[0], &adjncy[0],
						NULL, options, &perm[0], &iperm[0]);
					Q.resize(size);
					std::copy(perm.begin(), perm.end(), Q.begin());
				}
				else std::cout << "Unknown algorithm: " << mt << std::endl;
#else //USE_METIS
				std::cout << "No METIS library support" << std::endl;
#endif //USE_METIS
			}
			if (!part.empty())
			{
				//add another block for separator
				if (GetParameters().template Get<int>("separator"))
				{
#if 0
					for (idx_t i = 0; i < size; ++i)
					{
						for (idx_t jt = 0; jt < G.RowSize(i) && part[i] != nparts; ++jt)
						{
							idx_t j = G.Col(i, jt);
							if (part[j] != nparts && part[i] != part[j])
							{
								part[i] = nparts;
								break;
							}
						}
					}
#else
					PriorityQueue<idx_t, std::greater<idx_t> > queue(G.Size());
					for (idx_t i = 0; i < size; ++i)
					{
						idx_t wgt = 0;
						for (idx_t jt = 0; jt < G.RowSize(i); ++jt)
							if (part[i] != part[G.Col(i, jt)])
								wgt++;
						if (wgt) queue.Push(i, wgt);
					}
					while (!queue.Empty())
					{
						idx_t i = queue.Pop();
						idx_t wgt = queue.GetKey(i);
						if (wgt)
						{
							for (idx_t jt = 0; jt < G.RowSize(i); ++jt)
							{
								idx_t j = G.Col(i, jt);
								if (queue.Contains(j) && part[i] != part[j])
									queue.DecreaseKey(j, queue.GetKey(j) - 1);
							}
							part[i] = nparts;
						}
					}
#endif
					nparts++; // another part for separator
				}
				//compute ordering
				idx_t cnt = 0;
				Q.resize(size);
				blocks.push_back(0);
				for (idx_t k = 0; k < nparts; ++k)
				{
					for (idx_t i = 0; i < size; ++i) if (part[i] == k)
						Q[i] = cnt++;
					std::cout << "part " << k << " end " << cnt << " size " << cnt - blocks.back() << std::endl;
					blocks.push_back(cnt);
				}
			}
		}
		{ //inverse ordering and assemble matrix
			std::vector<idx_t> iQ(size);
			for (idx_t i = 0; i < size; ++i) iQ[Q[i]] = i;
			//assemble B
			for (idx_t iB = 0; iB < size; ++iB)
			{
				idx_t i = iQ[iB];
				for (idx_t j = 0; j < A.RowSize(i); ++j)
					B.PushBack(Q[A.Col(i, j)], A.Val(i, j));
				B.FinalizeRow();
			}
		}
		if( sorted )
		{
			if( print )
				std::cout << "Sort row entries" << std::endl;
			B.SortRows();
		}
		if( write_matrix )
		{
			if( level > 0 )
			{
				if (GetParameters().template Get<std::string>("write_format") == "bin")
				{
					if (print) std::cout << "save B" << level << ".bin" << std::endl;
					B.SaveBinary("B" + to_string(level) + ".bin");
				}
				else
				{
					if (print) std::cout << "save B" << level << ".mtx" << std::endl;
					B.Save("B" + to_string(level) + ".mtx");
				}
				if (print) std::cout << "save order" << level << ".txt" << std::endl;
				SaveVector("order" + to_string(level) + ".txt", Q);
				if (!blocks.empty())
				{
					if (print) std::cout << "save blocks" << level << ".txt" << std::endl;
					SaveVector("blocks" + to_string(level) + ".txt", blocks);
				}
				if (!part.empty())
				{
					if (print) std::cout << "save part" << level << ".txt" << std::endl;
					SaveVector("part" + to_string(level) + ".txt", part);
				}
			}
			else 
			{
				if (GetParameters().template Get<std::string>("write_format") == "bin")
				{
					if (print) std::cout << "save B.bin" << std::endl;
					B.SaveBinary("B.bin");
				}
				else
				{
					if (print) std::cout << "save B.mtx" << std::endl;
					B.Save("B.mtx");
				}
				if (print) std::cout << "save order.txt" << std::endl;
				SaveVector("order.txt", Q);
				if (!blocks.empty())
				{
					if (print) std::cout << "save blocks.txt" << std::endl;
					SaveVector("blocks.txt", blocks);
				}
				if (!part.empty())
				{
					if (print) std::cout << "save part.txt" << std::endl;
					SaveVector("part.txt", part);
				}
			}
		}
		S.SetParameters(GetParameters().SubParameters("Solver"));
		Qb.resize(size);
		y.resize(size);
		return S.Setup(B);
	}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const
	{
		idx_t size = B.Size();
		for (idx_t k = 0; k < size; ++k)
		{
			Qb[Q[k]] = b[k];
			y[Q[k]] = x[k];
		}
		bool ret = S.Solve(Qb, y);
		for (idx_t k = 0; k < size; ++k)
			x[k] = y[Q[k]];
		return ret;
	}
	size_t Bytes() const {return S.Bytes() + B.Bytes() + get_bytes(Qb) + get_bytes(y) + get_bytes(Q);}
	//for external use
	const std::vector<idx_t>& GetOrder() { return Q; }
};
#endif //_METIS_ORDERING_H
