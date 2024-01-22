#include "bicgstab.h"
#include "amg_ruge_stuben.h"
#include "gauss_seidel.h"
#include "ilduc.h"
#include "fixed_stress.h"

int main(int argc, char ** argv)
{
	typedef AMGRugeStuben<GaussSeidel, BICGSTAB<ILDUC> > AMG;
	BICGSTAB< FixedStress<AMG, ILDUC> > Solver;

	if (argc < 3)
	{
		std::cout << argv[0] << " N matrix.mtx [rhs.mtx] [sol.mtx]" << std::endl;
		Solver.GetParameters().Save("params_default.txt");
		Solver.GetParameters().SaveRaw("params_default.raw");
	}
	else
	{
		CSRMatrix A;
		std::vector<double> x,b;
		int N = atoi(argv[1]);
		if( std::string(argv[2]).find(".bin") != std::string::npos )
			A.LoadBinary(std::string(argv[2]));
		else
			A.Load(std::string(argv[2]));
		if( argc > 3 )
		{
			if( std::string(argv[3]).find(".bin") != std::string::npos )
				LoadVectorBinary(std::string(argv[3]),b);
			else
				LoadVector(std::string(argv[3]),b);
		}
		else b.resize(A.Size(),1.0);
		if( argc > 4 ) 
		{
			if( std::string(argv[4]).find(".bin") != std::string::npos )
				LoadVectorBinary(std::string(argv[4]),x);
			else
				LoadVector(std::string(argv[4]),x);
		}
		
		
		
		Solver.GetParameters().Load("params.txt");
		std::cout << "Loaded parameters: " << std::endl;
		std::cout << "Set input block size: " << N << std::endl;
		Solver.GetParameters().Set("Preconditioner:block_beg", N);
		Solver.GetParameters().Set("Preconditioner:block_end", A.Size());
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
