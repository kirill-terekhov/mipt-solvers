#include "csrmatrix.h"

int main(int argc, char ** argv)
{
	if( argc < 2 )
			std::cout << argv[0] << " matrix.mtx [tol=1.0e-5]" << std::endl;
	else
	{
		double tol = 1.0e-5;
		CSRMatrix A;
		A.Load(std::string(argv[1]));
		if( argc > 2 ) 
			tol = atof(argv[2]);
		if( A.Symmetric(tol) ) 
		{
			std::cout << "matrix is symmetric" << std::endl;
			return 0;
		}
		else std::cout << "matrix is not symmetric" << std::endl;
	}
	return 1;

}
