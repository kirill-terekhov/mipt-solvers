#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

int main(int argc, char** argv)
{
	if (argc < 2)
		std::cout << "Usage: " << argv[0] << " N" << std::endl;
	else
	{
		// N x N regular grid
		//   *------------*
		//   |            |
		//   |  p(i,j)    |
		//   |            |
		//   *------------*
		// p_xx: - p(i-1,j) + 2 p(i,j) - p(i+1,j)
		// p_xx + p_yy
		//            -1 p(i,j+1)
		//-1 p(i-1,j) +4 p(i,j)   - 1 p(i+1,j)
		//            -1 p(i,j-1)
		const double stencil[9] =
		{
			 0.0,-1.0, 0.0,
			-1.0, 0.0,-1.0,
			 0.0,-1.0, 0.0,
		};
#define Is(ii,jj) (((ii)+1)*3 + ((jj)+1))
		// Unknowns:
		// Total unknowns: 
		// p - N*M
		// All boundaries are zero dirichlet
		int N = atoi(argv[1]);
		int Np = N * N;
		double h = 1.0 / (double)N, h2 = pow(h, 2);
		std::vector<int> ia(1, 0), ja;
		std::vector<double> a, b(Np, h2);
#define Ip(i,j) ((i)*N + (j))
		//Fill biharmonic system for pressure
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
			{
				size_t ijpos = a.size();
				size_t bpos = ia.size() - 1;
				//reserve my entry
				ja.push_back(Ip(i, j));
				a.push_back(4.0);
				for(int ii = -1; ii <= 1; ++ii)
					for (int jj = -1; jj <= 1; ++jj) if( stencil[Is(ii,jj)] )
					{
						if (i + ii >= 0 && i + ii < N && j + jj >= 0 && j + jj < N)
						{
							ja.push_back(Ip(i + ii, j + jj));
							a.push_back(stencil[Is(ii, jj)]);
						}
					}
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