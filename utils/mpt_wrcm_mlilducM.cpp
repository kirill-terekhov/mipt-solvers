#include "bicgstab.h"
#include "mlilducM.h"
#include "maximal_transversal.h"
#include "wrcm.h"


template<typename T>
class MPTWRCM : public Methods
{
	MaximalTransversal< WeightedReverseCuthillMckee< T > > method;
public:
	static Parameters DefaultParameters()
	{
		Parameters ret;
		ret.Set("name", "MPTWRCM");
		ret.SubParameters("MPT") = MaximalTransversal<DummySolver>::DefaultParameters();
		ret.SubParameters("WRCM") = WeightedReverseCuthillMckee<DummySolver>::DefaultParameters();
		ret.SubParameters("Solver") = T::DefaultParameters();
		return ret;
	}
	MPTWRCM() { GetParameters() = DefaultParameters(); }
	bool Setup(const CSRMatrix& A)
	{
		method.GetParameters() = GetParameters().SubParameters("MPT");
		method.GetParameters().SubParameters("Solver") = GetParameters().SubParameters("WRCM");
		method.GetParameters().SubParameters("Solver").SubParameters("Solver") = GetParameters().SubParameters("Solver");
		return method.Setup(A);
	}
	bool Solve(const std::vector<double>& b, std::vector<double>& x) const
	{
		return method.Solve(b, x);
	}
	size_t Bytes() const { return method.Bytes(); }
};


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
		
		BICGSTAB< MPTWRCM< MLILDUC< MPTWRCM > > > Solver;
		
		Solver.GetParameters().Save("params_default.txt");
		Solver.GetParameters().SaveRaw("params_default.raw");
		std::cout << "Loading params_mlilduc.txt" << std::endl;
		Solver.GetParameters().Load("params_mlilduc.txt");
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
