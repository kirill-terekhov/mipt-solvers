#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char** argv)
{
	if (argc < 4)
		std::cout << "Usage: " << argv[0] << " x0.txt x.txt N [M = N]" << std::endl;
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
		std::string file0(argv[1]);
		std::string file(argv[2]);
		int N = atoi(argv[3]);
		int M = argc > 4 ? atoi(argv[4]) : N;
		int Np = N * M, Ns = N * M;
		//shift in matrix
		int Sp = 0, Ss = Np;
		double hx = 1.0 / (double)N;
		double hy = 1.0 / (double)M;
		std::vector<double> x0(Np + Ns, 0.0), x(Np + Ns, 0.0);
#define Ip(i,j) (Sp + (i)*M + (j))
#define Is(i,j) (Ss + (i)*M + (j))
		std::ifstream inp0(file0.c_str()), inp(file.c_str());
		int Nt0, Nt;
		inp0 >> Nt0;
		inp >> Nt;
		if (Nt != Np + Ns || Nt0 != Np + Ns)
		{
			std::cout << "incorrect size! expect " << Np + Ns;
			if (Nt0 != Np + Ns) std::cout << " in " << file0 << " found " << Nt0;
			if (Nt != Np + Ns) std::cout << " in " << file << " found " << Nt;
			std::cout << std::endl;
			return -1;
		}
		for (int k = 0; k < Nt; ++k)
		{
			inp0 >> x0[k];
			inp >> x[k];
			//std::cout << "k: " << x0[k] << " " << x[k] << " " << x0[k] - x[k] << std::endl;
		}
		inp.close();
		std::ofstream outp("grid.vtk");
		outp << "# vtk DataFile Version 2.0" << std::endl;
		outp << "Results" << std::endl;
		outp << "ASCII" << std::endl;
		outp << "DATASET STRUCTURED_POINTS" << std::endl;
		outp << "DIMENSIONS " << N << " " << N << " " << 1 << std::endl;
		outp << "ORIGIN 0 0 0" << std::endl;
		outp << "SPACING " << hx << " " << hx << " " << 1 << std::endl;
		outp << "POINT_DATA " << N * N << std::endl;
		outp << "SCALARS PRESSURE double" << std::endl;
		outp << "LOOKUP_TABLE default" << std::endl;
		for (int j = 0; j < N; ++j)
			for (int i = 0; i < N; ++i)
				outp << x0[Ip(i, j)] - x[Ip(i, j)] << std::endl;
		outp << "SCALARS OIL double" << std::endl;
		outp << "LOOKUP_TABLE default" << std::endl;
		for (int j = 0; j < N; ++j)
			for (int i = 0; i < N; ++i)
				outp << x0[Is(i, j)] - x[Is(i, j)] << std::endl;
		outp << "SCALARS WATER double" << std::endl;
		outp << "LOOKUP_TABLE default" << std::endl;
		for (int j = 0; j < N; ++j)
			for (int i = 0; i < N; ++i)
				outp << 1 - (x0[Is(i, j)] - x[Is(i, j)]) << std::endl;
		outp.close();
	}
	return 0;
}