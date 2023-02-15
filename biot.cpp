#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

int main(int argc, char** argv)
{
	if (argc < 2)
		std::cout << "Usage: " << argv[0] << " N [M = N] [mu = 1] [lambda = 1] [alpha = 1] [zeta = 1] [kappa = 1] [dt = 1] [rw = 1e-4]" << std::endl;
	else
	{
		// Biot system for poroelasticity:
		// - div( mu * grad(u) + mu * grad(u)^T + lambda * div(u)) + alpha * grad(p) = 0,
		// zeta * dt(p) + alpha * dt(div(u)) - div( kappa/nu * grad(p) ) = q.
		// 
		// reformulation:
		//  - mu * laplace(u) - (mu + lambda) * grad(div(u)) + alpha * grad(p) = 0,
		// zeta * dt(p) + alpha * dt(div(u)) - kappa/nu * laplacian(p) = q.
		// 
		// parameters:
		// mu - first Lame parameter (elasticity)
		// lambda - second Lame parameter (elasticity)
		// alpha - Biot coefficient
		// zeta - specific storage coefficient
		// kappa - media permeability (divided by viscosity)
		// dt - time step size
		// rw - well radius
		// 
		// N x M regular grid
		// N+1 - nodes in Ox direction
		// M+1 - nodes in Oy direction
		//   *---v(i,j+1)---*
		//   |              |
		// u(i,j) p(i,j) u(i+1,j)
		//   |              |
		//   *----v(i,j)----*
		// all normals at u:
		// ->
		// all normals at v:
		// ^
		// |
		// (i,j) - cell
		// i in [0,N-1]
		// j in [0,M-1]
		// Unknowns:
		// u - vertical faces, displacement in Ox direction
		// v - horizontal faces, displacement in Oy direction
		// p - pressure at cell center
		// Total unknowns: 
		// u - (N+1)*M
		// v - N*(M+1)
		// p - N*M
		// Boundary conditions:
		// fixed displacement at all boundaries
		// no-flow at all boundaries
		// Initial state:
		// zero displacement
		// zero pressure
		int N = atoi(argv[1]);
		int M = argc > 2 ? atoi(argv[2]): N;
		double mu = argc > 3 ? atof(argv[3]) : 1.0;
		double lambda = argc > 4 ? atof(argv[4]) : 1.0;
		double alpha = argc > 5 ? atof(argv[5]) : 1.0;
		double zeta = argc > 6 ? atof(argv[6]) : 1.0;
		double kappa = argc > 7 ? atof(argv[7]) : 1.0;
		double dt = argc > 8 ? atof(argv[8]) : 1.0;
		double rw = argc > 9 ? atof(argv[9]) : 1e-4;
		int Nu = (N+1) * M;
		int Nv = N * (M+1);
		int Np = N * M;
		//shift in matrix
		int Su = 0, Sv = Nu, Sp = Nu + Nv;
		double hx = 1.0 / (double)N, hx2 = hx * hx;
		double hy = 1.0 / (double)M, hy2 = hy * hy;
		double hxy = hx * hy;
		std::vector<int> ia(1, 0), ja;
		std::vector<double> a, b(Nu + Nv + Np, 0.0);
#define Iu(i,j) (Su + (i)*M + (j))
#define Iv(i,j) (Sv + (i)*(M+1) + (j))
#define Ip(i,j) (Sp + (i)*M + (j))
		//Fill laplacian and pressure grad for u:
		for(int i = 0; i <= N; ++i)
			for (int j = 0; j < M; ++j)
			{
				size_t ijpos = a.size();
				size_t bpos = ia.size() - 1;
				//reserve my entry
				ja.push_back(Iu(i,j));
				a.push_back(0.0);
				b[bpos] = 0.0;
				if (i == 0 || i == N)
				{
					// u(i,j) = 0
					a[ijpos] = 1.0; //left side - Dirichlet
					b[bpos] = 0.0; //flow to the right
				}
				else
				{
					//left: laplacian plus u-part of divergence difference
					// (2 * mu + lambda) * (u(i,j) - u(i-1,j)) / hx2
					ja.push_back(Iu(i - 1, j));
					a.push_back(-(2 * mu + lambda) / hx2);
					a[ijpos] += (2 * mu + lambda) / hx2;
					//right: laplacian plus u-part of divergence difference
					// (2 * mu + lambda) * (u(i,j) - u(i+1,j)) / hx2
					ja.push_back(Iu(i + 1, j));
					a.push_back(-(2 * mu + lambda) / hx2);
					a[ijpos] += (2 * mu + lambda) / hx2;
					//bottom: laplacian only
					if (j == 0) //slip condition at bottom: mu * (u(i,j) - 0.0) / (0.5 * hy2)
						a[ijpos] += mu / (0.5 * hy2);
					else 
					{
						// mu * (u(i,j) - u(i,j-1)) / hy
						ja.push_back(Iu(i, j - 1));
						a.push_back(-mu / hy2);
						a[ijpos] += mu / hy2;
					}
					//top: laplacian only
					if (j == M - 1) //slip condition at top: mu * (u(i,j) - 0.0) / (0.5 * hy2)
						a[ijpos] += mu / (0.5 * hy2);
					else
					{
						// mu * (u(i,j) - u(i,j+1)) / hy
						ja.push_back(Iu(i, j + 1));
						a.push_back(-mu / hy2);
						a[ijpos] += mu / hy2;
					}
					// v-part of divergence difference
					// (mu + lambda) * (v(i,j) + v(i-1,j+1) - v(i,j+1) - v(i-1,j)) / (hx hy)
					ja.push_back(Iv(i, j));
					a.push_back((mu + lambda) / hxy);
					ja.push_back(Iv(i, j + 1));
					a.push_back(-(mu + lambda) / hxy);
					ja.push_back(Iv(i - 1, j));
					a.push_back(-(mu + lambda) / hxy);
					ja.push_back(Iv(i - 1, j + 1));
					a.push_back((mu + lambda) / hxy);
					// alpha * (p(i,j) - p(i-1,j))/hx
					//pressure from the left
					ja.push_back(Ip(i - 1, j));
					a.push_back(-alpha / hx);
					//pressure from the right
					ja.push_back(Ip(i, j));
					a.push_back(alpha / hx);
				}
				ia.push_back((int)ja.size()); //close row
			}
		//Fill laplacian and pressure grad for v:
		for (int i = 0; i < N; ++i)
			for (int j = 0; j <= M; ++j)
			{
				size_t ijpos = a.size();
				size_t bpos = ia.size() - 1;
				//reserve my entry
				ja.push_back(Iv(i, j));
				a.push_back(0.0);
				b[bpos] = 0.0;
				if (j == 0 || j == M) // slip condition
				{
					// v(i,j) = 0
					a[ijpos] = 1.0;
					b[bpos] = 0.0;
				}
				else
				{
					//bottom: laplacian plus v-part of divergence difference
					// (2 * mu + lambda) * (v(i,j) - v(i, j-1)) / hy
					ja.push_back(Iv(i, j - 1));
					a.push_back(-(2 * mu + lambda) / hy2);
					a[ijpos] += (2 * mu + lambda) / hy2;
					//top: laplacian plus v-part of divergence difference
					// (2 * mu + lambda) * (v(i,j) - v(i, j+1)) / hy
					ja.push_back(Iv(i, j + 1));
					a.push_back(-(2 * mu + lambda) / hy2);
					a[ijpos] += (2 * mu + lambda) / hy2;
					//left: laplacian only
					if (i == 0) //zero vertical at left: mu * (v(i,j) - 0.0) / (0.5*hx)
						a[ijpos] += mu / (0.5 * hx2);
					else
					{
						// mu * (v(i,j) - v(i-1,j)) / hx
						ja.push_back(Iv(i - 1, j));
						a.push_back(-mu / hx2);
						a[ijpos] += mu / hx2;
					}
					//right: laplacian only
					if(i == N-1) //zero vertical at right: mu * (v(i,j) - 0.0) / (0.5*hx)
						a[ijpos] += mu / (0.5 * hx2);
					else
					{
						// mu * (v(i,j) - v(i+1,j)) / hx
						ja.push_back(Iv(i + 1, j));
						a.push_back(-mu / hx2);
						a[ijpos] += mu / hx2;
					}
					// u-part of divergence difference
					// (mu + lambda) * (u(i,j) + u(i-1,j+1) - u(i,j+1) - u(i-1,j)) / (hx hy)
					ja.push_back(Iu(i, j));
					a.push_back((mu + lambda) / hxy);
					ja.push_back(Iu(i, j - 1));
					a.push_back(-(mu + lambda) / hxy);
					ja.push_back(Iu(i + 1, j));
					a.push_back(-(mu + lambda) / hxy);
					ja.push_back(Iu(i + 1, j - 1));
					a.push_back((mu + lambda) / hxy);
					// (p(i,j) - p(i,j-1)) / hy
					//pressure from the top
					ja.push_back(Ip(i, j));
					a.push_back(alpha / hy);
					//pressure from the bottom
					ja.push_back(Ip(i, j - 1));
					a.push_back(-alpha / hy);
				}
				ia.push_back((int)ja.size()); //close row
			}
		//Fill filtration equation for pressure
		double WI = 2 * M_PI * kappa / log(0.14 / rw * sqrt(hx2 + hy2));
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
			{
				size_t ijpos = a.size();
				size_t bpos = ia.size() - 1;
				//reserve my entry
				ja.push_back(Ip(i, j));
				a.push_back(0.0);
				b[bpos] = 0.0;
				if (i == N / 2 && j == M / 2) //add injection well with bottom hole pressure 100
				{
					// WI * (p - 100.0)
					a[ijpos] += WI * dt;
					b[bpos] += WI * dt * 100.0;
					std::cout << "well at " << i << " " << j << std::endl;
				}
				//add specific storage term: 
				// zeta * (p - p0) / dt, p0 = 0
				a[ijpos] += zeta;
				//laplacian of p:
				//left:
				if (i > 0) //left BC is no-flow
				{
					// kappa * (p(i,j) - p(i-1,j)) / hx2
					ja.push_back(Ip(i - 1, j));
					a.push_back(-kappa * dt / hx2);
					a[ijpos] += kappa * dt / hx2;
				}
				//right:
				if (i < N - 1) //right BC is no-flow
				{
					// kappa * (p(i,j) - p(i+1,j)) / hx2
					ja.push_back(Ip(i + 1, j));
					a.push_back(-kappa * dt / hx2);
					a[ijpos] += kappa * dt / hx2;
				}
				//bottom:
				if (j > 0) //bottom BC is no-flow
				{
					// kappa * (p(i,j) - p(i,j-1)) / hy2
					ja.push_back(Ip(i, j - 1));
					a.push_back(-kappa * dt / hy2);
					a[ijpos] += kappa * dt / hy2;
				}
				//top:
				if (j < M - 1) //top BC is no-flow
				{
					// kappa * (p(i,j) - p(i,j+1)) / hy2
					ja.push_back(Ip(i, j + 1));
					a.push_back(-kappa * dt / hy2);
					a[ijpos] += kappa * dt / hy2;
				}
				//u-derivative
				ja.push_back(Iu(i, j));
				a.push_back(-alpha / hx);
				ja.push_back(Iu(i + 1, j));
				a.push_back(alpha / hx);
				//v-derivative
				ja.push_back(Iv(i, j));
				a.push_back(-alpha / hy);
				ja.push_back(Iv(i, j + 1));
				a.push_back(alpha / hy);

				ia.push_back((int)ja.size()); //close row
			}
		//write matrix
		std::ofstream output("A.mtx"), vec("b.txt");
		output << "%%MatrixMarket matrix coordinate real general" << std::endl;
		output << ia.size()-1 << " " << ia.size()-1 << " " << ja.size() << std::endl;
		vec << ia.size() - 1 << std::endl;
		//output << std::scientific;
		//output.precision(16);
		for (size_t i = 0; i < ia.size() - 1; ++i)
		{
			for (int j = ia[i]; j < ia[i + 1]; ++j)
				output << i + 1 << " " << ja[j] + 1 << " " << a[j] << std::endl;
			vec << b[i] << std::endl;
		}
		output.close();
	}
	return 0;
}