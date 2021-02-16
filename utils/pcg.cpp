#include "pcg.h"
#include "chebyshev.h"

int main(int argc, char ** argv)
{

	PCG< Chebyshev > Solver;

	if (argc < 2)
	{
		std::cout << argv[0] << " matrix.mtx [rhs.mtx] [sol.mtx]" << std::endl;
		Solver.GetParameters().Save("params_default.txt");
		Solver.GetParameters().SaveRaw("params_default.raw");
	}
	else
	{
		CSRMatrix A;
		std::vector<double> x,b;
		if( std::string(argv[1]).find(".bin") != std::string::npos )
			A.LoadBinary(std::string(argv[1]));
		else
			A.Load(std::string(argv[1]));
		if( argc > 2 )
		{
			if( std::string(argv[2]).find(".bin") != std::string::npos )
				LoadVectorBinary(std::string(argv[2]),b);
			else
				LoadVector(std::string(argv[2]),b);
		}
		else b.resize(A.Size(),1.0);
		if( argc > 3 ) 
		{
			if( std::string(argv[3]).find(".bin") != std::string::npos )
				LoadVectorBinary(std::string(argv[3]),x);
			else
				LoadVector(std::string(argv[3]),x);
		}
		
		
		
		Solver.GetParameters().Load("params_cg.txt");
		std::cout << "Loaded parameters: " << std::endl;
		Solver.GetParameters().Print();
		if( Solver.Setup(A) && Solver.Solve(b,x) )
		{
			std::cout << "Final residual " << Resid(A,b,x) << std::endl;
			SaveVector(std::string("solution"),x);
			return 0;
		}
		else std::cout << "Solution failed" << std::endl;
	}
	return 1;
}
