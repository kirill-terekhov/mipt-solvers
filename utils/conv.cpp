#include "csrmatrix.h"




int main(int argc, char ** argv)
{
	if( argc < 2 )
		std::cout << argv[0] << " matrix.mtx [rhs.mtx] [sol.mtx]" << std::endl;
	else
	{
		CSRMatrix A;
		std::vector<double> x,b;
		std::cout << "Loading " << argv[1] << std::endl;
		if( std::string(argv[1]).find(".bin") != std::string::npos )
			A.LoadBinary(std::string(argv[1]));
		else
			A.Load(std::string(argv[1]));
		std::cout << "Loaded " << argv[1] << ", saving A_out.bin" << std::endl;
		A.SaveBinary("A_out.bin");
		std::cout << "A_out.bin saved" << std::endl;
		if( argc > 2 )
		{
			std::cout << "Loading " << argv[2] << std::endl;
			if( std::string(argv[2]).find(".bin") != std::string::npos )
				LoadVectorBinary(std::string(argv[2]),b);
			else
				LoadVector(std::string(argv[2]),b);
			std::cout << "Loaded " << argv[2] << ", saving b_out.bin" << std::endl;
			SaveVectorBinary("b_out.bin",b);
			std::cout << "b_out.bin saved" << std::endl;
		}
		if( argc > 3 ) 
		{
			std::cout << "Loading " << argv[3] << std::endl;
			if( std::string(argv[3]).find(".bin") != std::string::npos )
				LoadVectorBinary(std::string(argv[3]),x);
			else
				LoadVector(std::string(argv[3]),x);
			std::cout << "Loaded " << argv[3] << ", saving x_out.bin" << std::endl;
			SaveVectorBinary("x_out.bin",x);
			std::cout << "x_out.bin saved" << std::endl;
		}
		return 0;
	}
	return 1;
}
