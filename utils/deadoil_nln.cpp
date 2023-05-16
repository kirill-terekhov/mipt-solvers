#include <iostream>
#include <fstream>
#include <vector>
#include "amg_ruge_stuben.h"
#include "gauss_seidel.h"
#include "ilduc.h"
#include "cpr.h"
#include "two_stage.h"
#include "bicgstab.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

const double natol = 1.0e-9;
const double nrtol = 1.0e-5;

int main(int argc, char** argv)
{
	if (argc < 2)
		std::cout << "Usage: " << argv[0] << " N [M = N] [kappa = 1] [phi = 0.5] [dt = 0.1] [rw = 1.0e-5]" << std::endl;
	else
	{
		// N x M regular grid
		// N - cells in Ox direction
		// M - cells in Oy direction
		//   *--------------*
		//   |              |
		//   | p(i,j),s(i,j)|
		//   |              |
		//   *--------------*
		// (i,j) - cell
		// i in [0,N-1]
		// j in [0,M-1]
		// Unknowns:
		// p - water pressure at cell center
		// s - oil saturation at cell center
		// Initial:
		// p = 100 [bar]
		// s = 0.75 [1]
		// Total unknowns: 
		// p - N*M
		// S - N*M
		// Paramters:
		// kappa - permeability
		// phi - porosity
		// dt - time step size
		// Equations:
		// water balance: dt(phi (1-s)) - div(krw(s) kappa grad(p)) = qw
		// oil balance:   dt(phi s) - div(kro(s) kappa grad(p)) = qo
		// phi - constant
		// kro(s) = s^2
		// krw(s) = (1-s)^2
		int N = atoi(argv[1]);
		int M = argc > 2 ? atoi(argv[2]): N;
		double kappa = argc > 3 ? atof(argv[3]) : 1.0;
		double phi = argc > 4 ? atof(argv[4]) : 0.5;
		double dt = argc > 5 ? atof(argv[5]) : 0.1;
		double rw = argc > 6 ? atof(argv[6]) : 1.0e-5;
		int Np = N * M, Ns = N * M;
		//shift in matrix
		int Sp = 0, Ss = Np;
		double hx = 1.0 / (double)N, hx2 = hx * hx;
		double hy = 1.0 / (double)M, hy2 = hy * hy;
		double u, skr, kr, h;
		double mu[2] = { 0.01, 1 };
		std::vector<idx_t> ia(1, 0), ja;
		std::vector<double> a, b(Np + Ns, 0.0);
		std::vector<double> x(Np + Ns, 0.0), x0(Np + Ns, 0.0);
#define Ip(i,j) (Sp + (i)*M + (j))
#define Is(i,j) (Ss + (i)*M + (j))
#define p(i,j) (x[Ip((i),(j))])
#define s(i,j) (x[Is((i),(j))])
#define s0(i,j) (x0[Is((i),(j))])
		std::fill(x0.begin() + Sp, x0.begin() + Sp + Np, 100); //initial water pressure 100 bar
		std::fill(x0.begin() + Ss, x0.begin() + Ss + Ns, 0.75); //initial oil saturation 0.75
		std::copy(x0.begin(), x0.end(), x.begin()); //copy initial to current solution
		double WI = 2 * M_PI * kappa / log(0.14 / rw * sqrt(hx2 + hy2)); //well index
		int step = 0;
		const double T = 10.0; //total time
		for (double t = 0; t < T; t += dt)
		{
			std::cout << std::setw(10) << t << "|" << T << std::endl;
			double norm0 = 1;
			for (int k = 0; k < 20; ++k)
			{
				a.clear();
				ia.clear();
				ja.clear();
				ia.resize(1, 0);
				std::fill(b.begin(), b.end(), 0);
				//fill water-oil balance equation
				for (int q = 0; q < 2; ++q) // q = 0 - water, q = 1 - oil
					for (int i = 0; i < N; ++i)
						for (int j = 0; j < M; ++j)
						{
							size_t ijpos = a.size();
							size_t bpos = ia.size() - 1;
							//reserve my entry
							ja.push_back(Ip(i, j));
							a.push_back(0.0);
							//accumulation term -phi*dt(S): -phi*(s(i,j) - s0(i,j))
							size_t sijpos = a.size();
							ja.push_back(Is(i, j));
							a.push_back(phi * (2 * q - 1));
							b[bpos] += phi * (2 * q - 1) * (s(i, j) - s0(i, j));
							//transport term -div( kr(s) kappa grad(p) )
							for (int ii = -1; ii <= 1; ++ii)
								for (int jj = -1; jj <= 1; ++jj)
									if (abs(ii) + abs(jj) == 1) //only one of the directions
									{
										h = (abs(ii) * hx2 + abs(jj) * hy2) / dt;
										if (i + ii >= 0 && i + ii < N && j + jj >= 0 && j + jj < M) //skip no-flow boundaries
										{
											//velocity: kappa * (p(i,j) - p(i+ii,j+jj))/h
											u = kappa / mu[q] * (p(i, j) - p(i + ii, j + jj)) / h;
											if (u > 0) //upstream
											{
												skr = 1 - q + (2 * q - 1) * s(i, j); //water: 1-s, oil: s
												skr = std::min(std::max(skr, 0.0), 1.0);
												kr = pow(skr, 2);
												a[sijpos] += 2.0 * (2 * q - 1) * skr * u;
											}
											else //downstream
											{
												skr = 1 - q + (2 * q - 1) * s(i + ii, j + jj); //water: 1-s, oil: s
												skr = std::min(std::max(skr, 0.0), 1.0);
												kr = pow(skr, 2);
												if (u)
												{
													ja.push_back(Is(i + ii, j + jj));
													a.push_back(2.0 * (2 * q - 1) * skr * u);
												}
											}
											ja.push_back(Ip(i + ii, j + jj));
											a.push_back(-kr * kappa / mu[q] / h);
											a[ijpos] += kr * kappa / mu[q] / h;
											b[bpos] += kr * u;
										}
									}
							if (i == 0 && j == 0 && q == 0) //add water injection well with bottom hole pressure 150 [bhp]
							{
								// WI * (p - 10)
								u = WI * dt * (p(i, j) - 150);
								if (u > 0)
								{
									std::cout << "upstream in injector: " << u << std::endl;
								}
								kr = 1;
								a[ijpos] += kr * WI * dt;
								b[bpos] += kr * u;
								//std::cout << "well at " << i << " " << j << " WI " << WI << " kr " << kr << " add " << kr * WI * dt << std::endl;
							}
							if (i == N - 1 && j == M - 1) //add production well with bottom hole pressure 10 [bhp]
							{
								// WI * (p - 10)
								u = WI * dt * (p(i, j) - 10);
								if (u > 0) //upstream
								{
									skr = 1 - q + (2 * q - 1) * s(i, j); //water: 1-s, oil: s
									skr = std::min(std::max(skr, 0.0), 1.0);
									kr = pow(skr, 2);
									a[sijpos] += 2.0 * (2 * q - 1) * skr * u;
								}
								else if (u) //downstream
								{
									std::cout << "downstream in producer: " << u << std::endl;
									skr = 1 - q + (2 * q - 1) * s0(i, j); //water: 1-s, oil: s
									skr = std::min(std::max(skr, 0.0), 1.0);
									kr = pow(skr, 2);
								}
								a[ijpos] += kr * WI * dt;
								b[bpos] += kr * u;
								//std::cout << "well at " << i << " " << j << " WI " << WI << " kr " << kr << " add " << kr * WI * dt << std::endl;
							}
							ia.push_back((int)ja.size()); //close row
						}
				//check convergence
				double norm = 0.0;
				for (size_t l = 0; l < b.size(); ++l)
					norm += b[l] * b[l];
				norm = sqrt(norm);
				if (k == 0) norm0 = norm;
				if (norm < natol || norm < nrtol * norm0)
				{
					std::cout << "converged with " << norm << " initial " << norm0 << std::endl;
					break;
				}
				else std::cout << "norm " << norm << " relative " << norm / norm0;// << std::endl;
				//solve system
				{
					CSRMatrix A(ia, ja, a);
					//BICGSTAB< MaximalTransversal< WeightedReverseCuthillMckee< ILDUC > > > Solver;
					BICGSTAB< CPR< AMGRugeStuben<GaussSeidel, BICGSTAB<ILDUC> >, GaussSeidel, TwoStage> > Solver;
					std::vector<double> dx;
					Solver.GetParameters().SetRecursive("verbosity", "0");
					Solver.GetParameters().Set("Preconditioner:block_beg", 0);
					Solver.GetParameters().Set("Preconditioner:block_end", N* M);
					//Solver.GetParameters().Set("Preconditioner:Solver:Solver:inverse_estimation", "1");
					//Solver.GetParameters().Set("Preconditioner:Solver:Solver:drop_tolerance", "0.001");
					if (Solver.Setup(A) && Solver.Solve(b, dx))
					{
						double alpha = 1.0, ds;
						for(int i = 0; i < N; ++i)
							for (int j = 0; j < M; ++j)
							{
								ds = -dx[Is(i, j)];
								if (s(i, j) + ds < 0.0)
									alpha = std::min(alpha, (0.0 - s(i, j)) / ds);
								if (s(i, j) + ds > 1.0)
									alpha = std::min(alpha, (1.0 - s(i, j)) / ds);
							}
						std::cout << " update alpha: " << alpha << std::endl;
						for (size_t l = 0; l < x.size(); ++l)
							x[l] -= alpha * dx[l];
					}
					else
					{
						std::cout << "Matrix not solved!" << std::endl;
						return -1;
					}
				}
			}
			//advance time
			std::copy(x.begin(), x.end(), x0.begin()); //take current solution as next step initial solution

			//record grid
			{
				std::ofstream outp("grid" + std::to_string(step)+".vtk");
				outp << "# vtk DataFile Version 2.0" << std::endl;
				outp << "Results" << std::endl;
				outp << "ASCII" << std::endl;
				outp << "DATASET STRUCTURED_POINTS" << std::endl;
				outp << "DIMENSIONS " << N << " " << M << " " << 1 << std::endl;
				outp << "ORIGIN 0 0 0" << std::endl;
				outp << "SPACING " << hx << " " << hy << " " << 1 << std::endl;
				outp << "POINT_DATA " << N * M << std::endl;
				outp << "SCALARS PRESSURE double" << std::endl;
				outp << "LOOKUP_TABLE default" << std::endl;
				for (int j = 0; j < M; ++j)
					for (int i = 0; i < N; ++i)
						outp << p(i, j) << std::endl;
				outp << "SCALARS OIL double" << std::endl;
				outp << "LOOKUP_TABLE default" << std::endl;
				for (int j = 0; j < M; ++j)
					for (int i = 0; i < N; ++i)
						outp << s(i, j) << std::endl;
				outp << "SCALARS WATER double" << std::endl;
				outp << "LOOKUP_TABLE default" << std::endl;
				for (int j = 0; j < M; ++j)
					for (int i = 0; i < N; ++i)
						outp << 1.0 - s(i, j) << std::endl;
				outp.close();
				step++;
			}
		}
	}
	return 0;
}