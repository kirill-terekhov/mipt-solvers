#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char** argv)
{
	if (argc < 3)
		std::cout << "Usage: " << argv[0] << " x.txt N [M = N]" << std::endl;
	else
	{
		// N x N regular grid
		//   *---------*
		//   |         |
		//   |  p(i,j) |
		//   |         |
		//   *---------*
		// (i,j) - cell
		// i in [0,N-1]
		// j in [0,N-1]
		// Unknowns:
		// p - pressure at cell center
		// Total unknowns: 
		// p - N*N
		std::string file(argv[1]);
		int N = atoi(argv[2]);
		int Np = N * N;
		//shift in matrix
		double hx = 1.0 / (double)N;
		std::vector<double> x;
#define Ip(i,j) ((i)*N + (j))
		std::ifstream inp(file.c_str());
		int Nt;
		inp >> Nt;
		if (Nt != Np)
		{
			std::cout << "incorrect size!" << std::endl;
			return -1;
		}
		x.resize(Nt);
		for (int k = 0; k < Nt; ++k)
			inp >> x[k];
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
				outp << x[Ip(i, j)] << std::endl;
		outp.close();
	}
	return 0;
}