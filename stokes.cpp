#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char** argv)
{
	if (argc < 2)
		std::cout << "Usage: " << argv[0] << " N [M = N] [visc = 1]" << std::endl;
	else
	{
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
		// u - vertical faces, velocity in Ox direction
		// v - horizontal faces, velocity in Oy direction
		// p - pressure at cell center
		// Total unknowns: 
		// u - (N+1)*M
		// v - N*(M+1)
		// p - N*M
		// left side - inflow with velocity u = 1
		// right side - outflow with pressure p = 0
		int N = atoi(argv[1]);
		int M = argc > 2 ? atoi(argv[2]): N;
		double mu = argc > 3 ? atof(argv[3]) : 1.0;
		int Nu = (N+1) * M;
		int Nv = N * (M+1);
		int Np = N * M;
		//shift in matrix
		int Su = 0, Sv = Nu, Sp = Nu + Nv;
		double hx = 1.0 / (double)N, hx2 = hx * hx;
		double hy = 1.0 / (double)M, hy2 = hy * hy;
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
				if (i == 0)
				{
					// u(i,j) = 1
					a[ijpos] = 1.0; //left side - Dirichlet
					b[bpos] = 1.0; //flow to the right
				}
				else
				{
					//left
					// mu * (u(i,j) - u(i-1,j)) / hx
					ja.push_back(Iu(i - 1, j));
					a.push_back(-mu / hx2);
					a[ijpos] += mu / hx2;
					//bottom 
					if (j == 0) //slip condition at bottom: mu * (u(i,j) - 0.0) / (0.5 * hy)
						a[ijpos] += mu / (0.5 * hy2);
					else 
					{
						// mu * (u(i,j) - u(i,j-1)) / hy
						ja.push_back(Iu(i, j - 1));
						a.push_back(-mu / hy2);
						a[ijpos] += mu / hy2;
					}
					//top
					if (j == M - 1) //slip condition at top: mu * (u(i,j) - 0.0) / (0.5 * hy)
						a[ijpos] += mu / (0.5 * hy2);
					else
					{
						// mu * (u(i,j) - u(i,j+1)) / hy
						ja.push_back(Iu(i, j + 1));
						a.push_back(-mu / hy2);
						a[ijpos] += mu / hy2;
					}
					//right
					if (i < N) //outflow to the right at i == N
					{
						// mu * (u(i,j) - u(i+1,j)) / hx
						ja.push_back(Iu(i + 1, j));
						a.push_back(-mu / hx2);
						a[ijpos] += mu / hx2;
					}
					if (i == N)//outflow with zero p at right
					{
						// (0 - p(i-1,j)) / (0.5*hx)
						//pressure from the left
						ja.push_back(Ip(i - 1, j));
						a.push_back(-1.0 / (0.5 * hx));
					}
					else 
					{
						// (p(i,j) - p(i-1,j))/hx
						//pressure from the left
						ja.push_back(Ip(i - 1, j));
						a.push_back(-1.0 / hx);
						//pressure from the right
						ja.push_back(Ip(i, j));
						a.push_back(1.0 / hx);
					}
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
				if (j == 0 || j == M) // slip condition
				{
					// v(i,j) = 0
					a[ijpos] = 1.0;
					b[bpos] = 0.0;
				}
				else
				{
					//left
					if (i == 0) //zero vertical at left: mu * (v(i,j) - 0.0) / (0.5*hx)
						a[ijpos] += mu / (0.5 * hx2);
					else
					{
						// mu * (v(i,j) - v(i-1,j)) / hx
						ja.push_back(Iv(i - 1, j));
						a.push_back(-mu / hx2);
						a[ijpos] += mu / hx2;
					}
					//bottom
					// mu * (v(i,j) - v(i, j-1)) / hy
					ja.push_back(Iv(i, j - 1));
					a.push_back(-mu / hy2);
					a[ijpos] += mu / hy2;
					//top
					// mu * (v(i,j) - v(i, j+1)) / hy
					ja.push_back(Iv(i, j + 1));
					a.push_back(-mu / hy2);
					a[ijpos] += mu / hy2;
					//right
					if (i < N - 1) //outflow to the right at i == N
					{
						// mu * (v(i,j) - v(i+1,j)) / hx
						ja.push_back(Iv(i + 1, j));
						a.push_back(-mu / hx2);
						a[ijpos] += mu / hx2;
					}
					// (p(i,j) - p(i,j-1)) / hy
					//pressure from the top
					ja.push_back(Ip(i, j));
					a.push_back(1.0 / hy);
					//pressure from the bottom
					ja.push_back(Ip(i, j - 1));
					a.push_back(-1.0 / hy);
				}
				ia.push_back((int)ja.size()); //close row
			}
		//Fill divergence-free condition for pressure
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
			{
				size_t ijpos = a.size();
				size_t bpos = ia.size() - 1;
				//reserve my entry
				ja.push_back(Ip(i, j));
				a.push_back(0.0);
				b[bpos] = 0.0;
				//u-derivative
				ja.push_back(Iu(i, j));
				a.push_back(-1.0 / hx);
				ja.push_back(Iu(i + 1, j));
				a.push_back(1.0 / hx);
				//v-derivative
				ja.push_back(Iv(i, j));
				a.push_back(-1.0 / hy);
				ja.push_back(Iv(i, j + 1));
				a.push_back(1.0 / hy);

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