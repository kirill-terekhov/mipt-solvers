#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

const bool write_sol = false;
const double ndtol = 1.0e+20;
const double nrtol = 1.0e-13;
const double natol = 1.0e-10;
const double p = 2.0;
const double s = 1.0;
const int kmax = 1500;

int main(int argc, char** argv)
{
	if (argc < 1)
		std::cout << "Usage: " << argv[0] << " [N=50000] [epsilon1=0.3] [epsilon2=0.005] [depsilon=1.05]" << std::endl;
	else
	{
		// N x 1 regular grid
		// N+1 - nodes in Ox direction
		// equation: eps * laplace(p) + 2*(1-x^2)*p + p^2 = 1
		// left side - zero dirichlet: p(-1) = 0
		// right side - zero dirichlet: p(1) = 0
		size_t N = argc > 1 ? atoi(argv[1]) : 10000;
		double eps1 = argc > 2 ? atof(argv[2]) : 0.3;
		double eps2 = argc > 3 ? atof(argv[3]) : 0.005;
		double deps = argc > 4 ? atof(argv[4]) : 1.05;
		double h = 2.0 / (N + 1);
		std::cout << "Epsilon " << eps1 << ":" << eps2 << " dEpsilon " << deps << std::endl;
		std::vector<double> al(N), ad(N), ar(N), b(N, 0.0), x(N), dx(N);
		std::ofstream fout("out.txt");
		double eps = eps1;
		std::vector< std::vector<double> > sol0, sol1;
		std::vector<int> its;
		std::vector<bool> valid;
		sol0.push_back(std::vector<double>(N, 1.0));
		sol0.push_back(std::vector<double>(N, 0.0));
		int step = 0;
		while ((eps >= eps2 && eps1 > eps2) || (eps < eps2 && eps1 < eps2))
		{
			std::cout << "epsilon: " << eps << std::endl;
			double norm0 = 1;
			if (sol0.empty())
			{
				std::cout << "No initial solution!" << std::endl;
				break;
			}
			its.clear();
			valid.clear();
			valid.resize(sol0.size(), true);
			size_t nsols = sol0.size();
			while (nsols)
			{
				for (size_t nsol = 0; nsol < sol0.size(); ++nsol) if(valid[nsol])
				{
					std::copy(sol0[nsol].begin(), sol0[nsol].end(), x.begin());
					bool success = false;
					int it;
					for (int k = 0; k < kmax; ++k)
					{
						it = k;
						std::fill(b.begin(), b.end(), 0);
						std::fill(al.begin(), al.end(), 0);
						std::fill(ad.begin(), ad.end(), 0);
						std::fill(ar.begin(), ar.end(), 0);
						//left boundary
						al[0] = 0.0;
						ad[0] = -2 * eps / h;
						ar[0] = eps / h;
						b[0] = eps * (0.0 + x[1] - 2 * x[0]) / h;
						//right boundary
						al[N - 1] = eps / h;
						ad[N - 1] = -2 * eps / h;
						ar[N - 1] = 0.0;
						b[N - 1] = eps * (0.0 + x[N - 2] - 2 * x[N - 1]) / h;
						//internal
						for (size_t i = 1; i < N - 1; ++i)
						{
							al[i] = eps / h;
							ad[i] = -2 * eps / h;
							ar[i] = eps / h;
							b[i] = eps * (x[i - 1] + x[i + 1] - 2 * x[i]) / h;
						}
						//remaining part
						for (size_t i = 0; i < N; ++i)
						{
							// (N-1+1)*h = 2 - h
							// (N+1) * h = 2
							// h = 2/(N+1)
							double c = -1 + (i + 1) * h;
							b[i] += (2 * (1 - c * c) * x[i] + x[i] * x[i] - 1) * h;
							ad[i] += (2 * (1 - c * c) + 2 * x[i]) * h;
						}
						//check convergence
						double norm = 0.0;
						for (size_t l = 0; l < b.size(); ++l)
							norm += b[l] * b[l];
						norm = sqrt(norm / N);
						//deflation
						for (size_t ndef = 0; ndef < sol1.size(); ++ndef)
						{
							double M = 0.0;
							for (size_t k = 0; k < x.size(); ++k)
								M += (x[k] - sol1[ndef][k]) * (x[k] - sol1[ndef][k]);
							M = 1.0 / pow(M, p / 2.0) + s;
							norm *= M;
						}
						if (k == 0) norm0 = norm;
						if (norm < natol || norm < nrtol * norm0)
						{
							std::cout << "converged with " << std::setw(12) << norm 
								<< " initial " << std::setw(12) << norm0 
								<< " steps " << std::setw(4) << k
								<< " from " << std::setw(4) << nsol << std::endl;
							success = true;
							break;
						}
						if (norm > ndtol || k == kmax - 1)
						{
							std::cout << "diverged with " << std::setw(12) << norm
								<< " initial " << std::setw(12) << norm0
								<< " steps " << std::setw(4) << k
								<< " from " << std::setw(4) << nsol << std::endl;
							break;
						}
						//solve system
						{
							//tridiagonal system solver
							for (size_t i = 1; i < N; ++i)
							{
								double w = al[i] / ad[i - 1];
								ad[i] -= w * ar[i - 1];
								b[i]  -= w *  b[i - 1];
							}
							dx[N - 1] = b[N - 1] / ad[N - 1];
							for (size_t it = N - 1; it > 0; --it)
							{
								size_t i = it - 1;
								dx[i] = (b[i] - ar[i] * dx[i + 1]) / ad[i];
							}
							
							{
								//deflation
								double Edx0 = 0.0;
								for (size_t ndef = 0; ndef < sol1.size(); ++ndef)
								{
									double M = 0, Edx = 0;
									for (size_t k = 0; k < sol1[ndef].size(); ++k)
									{
										M += (x[k] - sol1[ndef][k]) * (x[k] - sol1[ndef][k]);
										Edx += p * (x[k] - sol1[ndef][k]) * dx[k];
									}
									Edx /= pow(M, p / 2.0 + 1.0);
									M = 1.0 / pow(M, p / 2.0) + s;
									Edx0 += Edx / M;
								}
								for (size_t k = 0; k < dx.size(); ++k)
									x[k] -= (1.0 + Edx0 / (1.0 - Edx0)) * dx[k];
							}
						}
					}
					if (success) //add solution
					{
						sol1.push_back(x);
						its.push_back(it);
					}
					else
					{
						valid[nsol] = false;
						nsols--;
					}
				}
			}
			for (size_t ndef = 0; ndef < sol1.size(); ++ndef)
			{
				double integr = 0.0;
				for (size_t k = 0; k < sol1[ndef].size(); ++k)
					integr += sol1[ndef][k] * sol1[ndef][k] * h;
				std::cout << ndef 
					<< " integral: " << std::setw(12) << integr 
					<< " derivative: " << std::setw(12) << (sol1[ndef][0] - 0.0) / h
					<< " output: " << std::setw(12) << integr * (sol1[ndef][0] - 0.0) / h
					<< " iterations: " << std::setw(4) << its[ndef] << std::endl;
				fout << std::setw(5) << eps << ";" << std::setw(12) << integr * (sol1[ndef][0] - 0.0) / h << std::endl;
			}
			if(write_sol)
			{
				std::ofstream sols("sols" + std::to_string(step) + ".txt");
				for (size_t k = 0; k < N; ++k)
				{
					sols << std::setw(12) << -1 + (k + 1) * h;
					for (size_t ndef = 0; ndef < sol1.size(); ++ndef)
						sols << ";" << std::setw(12) << sol1[ndef][k];
					sols << std::endl;
				}
				sols.close();
				std::cout << "written file sols" << step << ".txt" << std::endl;
			}
			std::swap(sol0, sol1);
			sol1.clear();
			eps /= deps;
			step++;
		}
	}
	return 0;
}
