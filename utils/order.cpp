#include "csrmatrix.h"
#include "metis_ordering.h"

int main(int argc, char ** argv)
{
	if( argc < 3 )
		std::cout << argv[0] << " nparts matrix.mtx [rhs.mtx] [sol.mtx]" << std::endl;
	else
	{
		idx_t nparts = atoi(argv[1]);
		std::vector<idx_t> Q;
		CSRMatrix C;
		{
			Metis<DummySolver> order;
			order.GetParameters().Set("parts", nparts);
			//order.GetParameters().Set("separator", "0");
			order.GetParameters().Set("write_matrix", "1");
			//order.GetParameters().Set("write_format", "bin");
			order.GetParameters().Set("debug", "0");
			CSRMatrix A;
			std::cout << "Loading matrix " << argv[2] << std::endl;
			if (std::string(argv[2]).find(".bin") != std::string::npos)
				A.LoadBinary(std::string(argv[2]));
			else
				A.Load(std::string(argv[2]));
			std::cout << "Loaded " << argv[2] << std::endl;
			order.Setup(A); //matrix is recorded inside
			Q = order.GetOrder(); //for vectors
			std::cout << "Done" << std::endl;
		}
		if( argc > 3 )
		{
			std::vector<double> b, Qb;
			std::cout << "Loading " << argv[3] << std::endl;
			if( std::string(argv[3]).find(".bin") != std::string::npos )
				LoadVectorBinary(std::string(argv[3]), b);
			else
				LoadVector(std::string(argv[3]), b);
			if (Q.size() != b.size())
				std::cout << "Error, vector " << argv[3] << " size " << b.size() << " mismatch ordering size " << Q.size() << std::endl;
			else
			{
				std::cout << "Loaded " << argv[3] << std::endl;
				Qb.resize(b.size());
				for (idx_t k = 0; k < b.size(); ++k)
					Qb[Q[k]] = b[k];
				SaveVectorBinary("Qb_out.bin",Qb);
				std::cout << "Qb_out.bin saved" << std::endl;
			}
		}
		if( argc > 4 ) 
		{
			std::vector<double> x, Qx;
			std::cout << "Loading " << argv[4] << std::endl;
			if( std::string(argv[4]).find(".bin") != std::string::npos )
				LoadVectorBinary(std::string(argv[4]), x);
			else
				LoadVector(std::string(argv[4]), x);
			if (Q.size() != x.size())
				std::cout << "Error, vector " << argv[4] << " size " << x.size() << " mismatch ordering size " << Q.size() << std::endl;
			else
			{
				std::cout << "Loaded " << argv[4] <<  std::endl;
				Qx.resize(x.size());
				for (idx_t k = 0; k < x.size(); ++k)
					Qx[Q[k]] = x[k];
				SaveVectorBinary("Qx_out.bin",Qx);
				std::cout << "x_out.bin saved" << std::endl;
			}
		}
		return 0;
	}
	return 1;
}
