#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char** argv)
{
	if (argc < 3)
		std::cout << "Usage: " << argv[0] << " x.txt N [M = N]" << std::endl;
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
		std::string file(argv[1]);
		int N = atoi(argv[2]);
		int M = argc > 3 ? atoi(argv[3]): N;
		int Nu = (N+1) * M;
		int Nv = N * (M+1);
		int Np = N * M;
		//shift in matrix
		int Su = 0, Sv = Nu, Sp = Nu + Nv;
		double hx = 1.0 / (double)N;
		double hy = 1.0 / (double)M;
		std::vector<double> x;
#define Iu(i,j) (Su + (i)*M + (j))
#define Iv(i,j) (Sv + (i)*(M+1) + (j))
#define Ip(i,j) (Sp + (i)*M + (j))
		std::ifstream inp(file.c_str());
		int Nt;
		inp >> Nt;
		if (Nt != Nu + Nv + Np)
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
		outp << "DIMENSIONS " << N << " " << M << " " << 1 << std::endl;
		outp << "ORIGIN 0 0 0" << std::endl;
		outp << "SPACING " << hx << " " << hy << " " << 1 << std::endl;
		outp << "POINT_DATA " << N * M << std::endl;
		outp << "SCALARS PRESSURE double" << std::endl;
		outp << "LOOKUP_TABLE default" << std::endl;
		for (int j = 0; j < M; ++j)
			for (int i = 0; i < N; ++i)
				outp << x[Ip(i, j)] << std::endl;
		outp << "SCALARS VELOCITY double 3" << std::endl;
		outp << "LOOKUP_TABLE default" << std::endl;
		for (int j = 0; j < M; ++j)
			for (int i = 0; i < N; ++i)
				outp << (x[Iu(i + 1, j)] + x[Iu(i, j)]) * 0.5 << " " << (x[Iv(i, j + 1)] + x[Iv(i, j)]) * 0.5 << " " << 0.0 << std::endl;
		outp << "SCALARS DIVERGENCE double" << std::endl;
		outp << "LOOKUP_TABLE default" << std::endl;
		for (int j = 0; j < M; ++j)
			for (int i = 0; i < N; ++i)
				outp << (x[Iu(i + 1, j)] - x[Iu(i, j)]) / hx + (x[Iv(i, j + 1)] - x[Iv(i, j)]) / hy << std::endl;
		outp.close();
	}
	return 0;
}