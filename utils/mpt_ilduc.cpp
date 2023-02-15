#include "bicgstab.h"
#include "ilduc.h"
#include "maximal_transversal.h"

int main(int argc, char ** argv)
{
	if( argc < 3 )
		std::cout << argv[0] << " matrix.mtx rhs.mtx [sol.mtx]" << std::endl;
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
		
		BICGSTAB< MaximalTransversal< ILDUC > > Solver;
		
		Solver.GetParameters().Save("params_default.txt");
		Solver.GetParameters().SaveRaw("params_default.raw");
		std::cout << "Loading params_mpt_ilduc.txt" << std::endl;
		Solver.GetParameters().Load("params_mpt_ilduc.txt");
		std::cout << "Loaded parameters: " << std::endl;
		Solver.GetParameters().Print();
		bool success;
		
		
		
		success = Solver.Setup(A);
		
		
		if( success )
			success = Solver.Solve(b,x);
		
		if( success )
		{
			std::cout << "Final residual " << Resid(A,b,x) << std::endl;
			SaveVector(std::string("solution"),x);
			return 0;
		}
		else std::cout << "Solution failed" << std::endl;
	}
	return 1;
}
