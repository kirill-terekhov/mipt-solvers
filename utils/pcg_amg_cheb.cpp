#include "pcg.h"
#include "amg_ruge_stuben.h"
#include "chebyshev.h"
#include "ilduc.h"

int main(int argc, char ** argv)
{

	PCG< AMGRugeStuben< Chebyshev, PCG<ILDUC> > > Solver;

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
		
		
		
		Solver.GetParameters().Load("params.txt");
		std::cout << "Loaded parameters: " << std::endl;
		Solver.GetParameters().Print();
		bool success = true;
		double t1, t2, t3, t4;
		t1 = get_time();
		success &= Solver.Setup(A);
		t2 = get_time();
		std::cout << "Setup time: " << t2 - t1 << std::endl;
		t3 = get_time();
		success &= Solver.Solve(b, x);
		t4 = get_time();
		std::cout << "Solve time: " << t2 - t1 << std::endl;
		if (success)
		{
			std::cout << "Time setup " << t2 - t1 << " iterations " << t4 - t3 << " solve " << t4 - t1 << std::endl;
			std::cout << "Solver consumed: " << Solver.Bytes() / 1024 << " KB" << std::endl;
			std::cout << "Matrix consumed: " << A.Bytes() / 1024 << " KB" << std::endl;
			std::cout << "Vector consumed: " << get_bytes(b) / 1024 << " KB" << std::endl;
			std::cout << "Final residual " << Resid(A,b,x) << std::endl;
			SaveVector(std::string("solution"),x);
			return 0;
		}
		else std::cout << "Solution failed" << std::endl;
	}
	return 1;
}
