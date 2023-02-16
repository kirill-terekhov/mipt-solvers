#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

int main(int argc, char** argv)
{
	if (argc < 2)
		std::cout << "Usage: " << argv[0] << " N [M = N] [kappa = 1] [phi = 0.5] [dt = 1] [rw = 1.0e-4]" << std::endl;
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
		double dt = argc > 5 ? atof(argv[5]) : 1.0;
		double rw = argc > 6 ? atof(argv[6]) : 1.0e-4;
		int Np = N * M, Ns = N * M;
		//shift in matrix
		int Sp = 0, Ss = Np;
		double hx = 1.0 / (double)N, hx2 = hx * hx;
		double hy = 1.0 / (double)M, hy2 = hy * hy;
		double u, skr, kr, h;
		std::vector<int> ia(1, 0), ja;
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
		//fill water-oil balance equation
		double WI = 2 * M_PI * kappa / log(0.14 / rw * sqrt(hx2 + hy2));
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
				for(int ii = -1; ii <= 1; ++ii)
					for(int jj = -1; jj <= 1; ++jj)
						if (abs(ii) + abs(jj) == 1) //only one of the directions
						{
							h = (abs(ii) * hx2 + abs(jj) * hy2) / dt;
							if (i + ii >= 0 && i + ii < N && j + jj >= 0 && j + jj < M) //skip no-flow boundaries
							{
								//velocity: kappa * (p(i,j) - p(i+ii,j+jj))/h
								u = kappa * (p(i, j) - p(i + ii, j + jj)) / h;
								if (u > 0) //upstream
								{
									skr = 1 - q + (2 * q - 1) * s(i, j); //water: 1-s, oil: s
									kr = pow(skr, 2);
									a[sijpos] += 2.0 * (2 * q - 1) * skr * u;
								}
								else //downstream
								{
									skr = 1 - q + (2 * q - 1) * s(i + ii, j + jj); //water: 1-s, oil: s
									kr = pow(skr, 2);
									if (u)
									{
										ja.push_back(Is(i + ii, j + jj));
										a.push_back(2.0 * (2 * q - 1) * skr * u);
									}
								}
								ja.push_back(Ip(i + ii, j + jj));
								a.push_back(-kr * kappa / h);
								a[ijpos] += kr * kappa / h;
								b[bpos] += kr * u;
							}
						}
				if (i == 0 && j == 0 && q == 0) //add water injection well with bottom hole pressure 150 [bhp]
				{
					// WI * (p - 10)
					u = WI * dt * (p(i, j) - 150);
					kr = 1;
					a[ijpos] += kr * WI * dt;
					b[bpos] += kr * u;
					std::cout << "well at " << i << " " << j << " WI " << WI << " kr " << kr << " add " << kr * WI * dt << std::endl;
				}
				if (i == N-1 && j == M-1) //add production well with bottom hole pressure 10 [bhp]
				{
					// WI * (p - 10)
					u = WI * dt * (p(i, j) - 10);
					if (u > 0) //upstream
					{
						skr = 1 - q + (2 * q - 1) * s(i, j); //water: 1-s, oil: s
						kr = pow(skr, 2);
						a[sijpos] += 2.0 * (2 * q - 1) * skr * u;
					}
					else if (u) //downstream
					{
						std::cout << "downstream in well: " << u << std::endl;
						skr = 1 - q + (2 * q - 1) * s0(i, j); //water: 1-s, oil: s
						kr = pow(skr, 2);
					}
					a[ijpos] += kr * WI * dt;
					b[bpos] += kr * u;
					std::cout << "well at " << i << " " << j << " WI " << WI << " kr " << kr << " add " << kr * WI * dt << std::endl;
				}
				ia.push_back((int)ja.size()); //close row
			}
		//write matrix
		std::ofstream output("A.mtx"), vec("b.txt"), vecx0("x0.txt");
		output << "%%MatrixMarket matrix coordinate real general" << std::endl;
		output << ia.size() - 1 << " " << ia.size() - 1 << " " << ja.size() << std::endl;
		vec << ia.size() - 1 << std::endl;
		vecx0 << ia.size() - 1 << std::endl;
		//output << std::scientific;
		//output.precision(16);
		for (size_t i = 0; i < ia.size() - 1; ++i)
		{
			for (int j = ia[i]; j < ia[i + 1]; ++j)
				output << i + 1 << " " << ja[j] + 1 << " " << a[j] << std::endl;
			vec << b[i] << std::endl;
			vecx0 << x0[i] << std::endl;
		}
		output.close();
	}
	return 0;
}