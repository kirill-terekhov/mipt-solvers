#ifndef _FIEDLER_H
#define _FIEDLER_H
#include "method.h"
#include "graph.h"
#include "priority_queue.h"
/*
 * Fiedler
 * 
 * compute fiedler vector from graph Laplacian 
 * using Housholder deflation and inverse power iteration
 * [1] An efficient and accurate method to compute the Fiedler
 * vector based on Householder deflation and inverse power
 * iteration
 * by Wu Jian-ping, Song Jun-qiang, Zhang Wei-min
 * TODO: AMG for inverse power iteration
 */
template<typename Solver>
class Fiedler : public Methods
{
	std::vector<idx_t> Q;
	CSRMatrix B;
	Solver S;
	mutable std::vector<double> Qb, y;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name","Fiedler");
		ret.Set("level","*");
		ret.Set("sorted", 1);
		ret.Set("verbosity",0); // 2 - output cg iterations
		ret.Set("separator", 1);
		ret.Set("parts", 4);
		ret.Set("maxiters", 50);
		ret.Set("tolerance", 1.0e-7);
		ret.Set("cg_iterations", 1000);
		ret.Set("cg_tolerance", 1.0e-9);
		ret.Set("write_matrix",1);
		ret.Set("write_format", "mtx");
		ret.SubParameters("Solver") = Solver::DefaultParameters();
		return ret;
	}
	Fiedler() {GetParameters() = DefaultParameters();}
	// Multiply y = (H L H) x
	//L - graph Laplacian
	//H - Housholder deflation represented by vectors u, v
	void HLHMult(const CSRMatrix& L, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& x, std::vector<double>& y) const
	{
		double a = Dot(u, x);
		double b = Dot(v, x);
		L.Multiply(1.0, x, 0.0, y);
		for (idx_t i = 0; i < L.Size(); ++i)
			y[i] -= u[i] * b + v[i] * a;
		y[0] = 0.0; // avoid error
	}
	//solve (H L H) x = b using conjugate gradient method
	//L - graph Laplacian
	//H - Housholder deflation represented by vectors u, v
	bool HLHSolve(const CSRMatrix& L, const std::vector<double> &u, const std::vector<double> &v, const std::vector<double>& b, std::vector<double>& x) const
	{
		idx_t N = L.Size();
		bool print = (GetParameters().template Get<int>("verbosity") & 2) ? true : false;
		int iters = 0, maxiters = GetParameters().template Get<int>("cg_iterations");
		double resid, beta, alpha, tol = GetParameters().template Get<double>("cg_tolerance");
		std::vector<double> r = b, p(x.size()), w(x.size());
		HLHMult(L, u, v, x, w);
		for (idx_t i = 0; i < N; ++i) r[i] -= w[i];
		p = r;
		resid = Dot(r, r);
		if (print)
			std::cout << "CG init " << std::setw(14) << sqrt(fabs(resid)) << " norm " << std::setw(14) << Norm(x) << std::endl;
		while (sqrt(fabs(resid)) > tol && iters < maxiters)
		{
			HLHMult(L, u, v, p, w);
			alpha = resid / Dot(w, p);
			for (idx_t i = 0; i < N; ++i) x[i] += alpha * p[i];
			for (idx_t i = 0; i < N; ++i) r[i] -= alpha * w[i];
			beta = 1.0 / resid;
			resid = Dot(r, r);
			beta *= resid;
			for (idx_t i = 0; i < N; ++i) p[i] = r[i] + beta * p[i];
			if (print)
				std::cout << "CG " << std::setw(4) << iters << " " << std::setw(14) << sqrt(fabs(resid)) << " norm " << std::setw(14) << Norm(x) << std::endl;
			iters++;
		}
		return sqrt(fabs(resid)) > tol ? false : true;
	}
	// Compute Fiedler vector for graph Laplacian L
	bool FindFiedler(const CSRMatrix& L, std::vector<double>& x, double & lambda) const
	{
		bool print = GetParameters().template Get<int>("verbosity") ? true : false;
		int maxiters = GetParameters().template Get<int>("maxiters");
		double tolerance = GetParameters().template Get<double>("tolerance");
		idx_t N = L.Size();
		std::vector<double> u(N), v(N), y(N);
		double alpha, beta, gamma, error, Lmax = 0;
		for(idx_t i = 0; i < N; ++i)	
			for(idx_t j = 0; j < L.RowSize(i); ++j)
				Lmax = std::max(L.Val(i,j), Lmax);
		// u = 1 + e0 / sqrt(N)
		std::fill(u.begin(), u.end(), 1.0);
		u[0] += sqrt(N);
		alpha = N + sqrt(N);
		// h = Lu / alpha
		L.Multiply(1.0, u, 0.0, v);
		for (idx_t i = 0; i < N; ++i)
			v[i] /= alpha;
		// gamma = u^T v / alpha
		gamma = Dot(u,v) / alpha;
		// v = h - gamma * u / 2
		for (idx_t i = 0; i < N; ++i)
			v[i] -= gamma * u[i] / 2.0;
		// set x = 1
		x.resize(N);
		std::fill(x.begin(), x.end(), 1.0);
		x[0] = 0.0;
		lambda = 1;
		for (int k = 0; k < maxiters; ++k)
		{
			beta = Dot(x, x); // beta = (x,x)
			for (idx_t i = 0; i < N; ++i) x[i] /= sqrt(beta);
			// y = (H L H) x
			HLHMult(L, u, v, x, y);
			//compute eigenvalue
			lambda = Dot(x, y);
			//compute error
			error = 0.0;
			for (idx_t i = 0; i < N; ++i)
				error = std::max(error, y[i] - lambda * x[i]);
			if( print ) 
				std::cout << "iter " << std::setw(3) << k << "|" << maxiters << " error " << std::setw(14) << error << " lambda " << std::setw(14) << lambda << std::endl;
			//check tolerance
			//if (error / Lmax < tolerance) break;
			if (error < tolerance) break;
			//solve for eigenvector
			// (H L H) y = x
			HLHSolve(L,u,v,x,y);
			x = y;
			//for (idx_t i = 0; i < N; ++i) std::cout << std::setw(14) << x[i];
			//std::cout << std::endl;
		}
		beta = Dot(u,x) / alpha;
		for (idx_t i = 0; i < N; ++i) x[i] -= beta * u[i];
		//return error / Lmax < tolerance ? true : false;
		return error < tolerance ? true : false;
	}
	//graph Laplacian of the matrix
	CSRMatrix GraphLaplacian(const CSRGraph & G) const
	{
		CSRMatrix L;
		for (idx_t i = 0; i < G.Size(); ++i)
		{
			bool diag = false;
			for (idx_t j = 0; j < G.RowSize(i); ++j)
			{
				if (G.Col(i, j) != i)
					L.PushBack(G.Col(i, j), -1.0);
				else diag = true;
			}
			L.PushBack(i, G.RowSize(i) - (diag ? 1.0 : 0.0));
			L.FinalizeRow();
		}
		return L;
	}
	void AddSeparator1(const CSRGraph& G, std::vector<int> & part, int & nparts)
	{
		for (idx_t i = 0; i < G.Size(); ++i)
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
		nparts++; // another part for separator
	}
	void AddSeparator(const CSRGraph& G, std::vector<int>& part, int& nparts)
	{
		PriorityQueue<idx_t, std::greater<idx_t> > queue(G.Size());
		for (idx_t i = 0; i < G.Size(); ++i)
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
		nparts++; // another part for separator
	}

	std::vector<idx_t> MakeParts(const std::vector<double>& x, double tolerance)
	{
		bool print = GetParameters().template Get<int>("verbosity") ? true : false;
		std::vector<idx_t> part(x.size());
		int nn = 0, nz = 0, np = 0;
		for (idx_t k = 0; k < (idx_t)x.size(); ++k)
		{
			if (x[k] > tolerance)
			{
				part[k] = 0;
				np++;
			}
			else if (x[k] < -tolerance)
			{
				part[k] = 1;
				nn++;
			}
			else
			{
				part[k] = rand() % 2;
				nz++;
			}
		}
		if (print)
			std::cout << "negative " << nn << " zero " << nz << " positive " << np << std::endl;
		return part;
	}

	CSRGraph GeneratePartGraph(const CSRMatrix& A, const std::vector<int> & part, int ipart)
	{
		CSRGraph G;
		std::vector<idx_t> id(A.Size());
		idx_t ind = 0;
		for (idx_t k = 0; k < A.Size(); ++k)
			if (part[k] == ipart) id[k] = ind++;
		for (idx_t k = 0; k < A.Size(); ++k) if( part[k] == ipart )
		{
			G.PushBack(id[k]);
			for (idx_t j = 0; j < A.RowSize(k); ++j)
				if (part[A.Col(k, j)] == ipart)
					G.PushBack(id[A.Col(k, j)]);
			G.FinalizeRow();
		}
		return G + G.Transpose();
	}

	bool Setup(const CSRMatrix& A)
	{
		bool print = GetParameters().template Get<int>("verbosity") ? true : false;
		bool write_matrix = GetParameters().template Get<int>("write_matrix") ? true : false;
		bool sorted = GetParameters().template Get<int>("sorted") ? true : false;
		int level = GetParameters().template Get<int>("level");
		int nparts = GetParameters().template Get<int>("parts");
		double tolerance = GetParameters().template Get<double>("tolerance");
		idx_t size = A.Size();
		//output data - which entity should belong to which part
		int depth = (int)ceil(log(nparts)/log(2.0));
		int nparts2 = pow(2, depth);
		if (print) std::cout << "nparts " << nparts << " log2 " << depth << " expected nparts " << nparts2 << std::endl;
		std::vector<int> part(A.Size(), 0);
		std::vector<idx_t> blocks;
		{ //memory
			double lambda = 0.0;
			std::vector<double> x;
			std::vector<idx_t> subpart(A.Size(), 0);
			CSRGraph G;
			CSRMatrix L;
			nparts = 1;
			for (int k = 0; k < depth; ++k)
			{
				for (int ipart = 0; ipart < pow(2, k); ++ipart)
				{
					//graph assembly
					G = GeneratePartGraph(A, part, ipart);
					L = GraphLaplacian(G);
					FindFiedler(L, x, lambda);
					subpart = MakeParts(x, tolerance * 2);
					idx_t ind = 0;
					for (idx_t k = 0; k < part.size(); ++k)
						if (part[k] == ipart)
						{
							if (subpart[ind] == 1)
								part[k] = nparts;
							ind++;
						}
					nparts++;
				}
				if (print)
					std::cout << "nparts: " << nparts << std::endl;
			}
			{
				//add another block for separator
				if (GetParameters().template Get<int>("separator"))
				{
					G = CSRGraph(A.get_ia(), A.get_ja());
					AddSeparator(G, part, nparts);
				}
				//compute ordering
				idx_t cnt = 0;
				Q.resize(size);
				blocks.push_back(0);
				for (idx_t k = 0; k < nparts; ++k)
				{
					for (idx_t i = 0; i < size; ++i) if (part[i] == k)
						Q[i] = cnt++;
					if (print) std::cout << "part " << k << " end " << cnt << " size " << cnt - blocks.back() << std::endl;
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
