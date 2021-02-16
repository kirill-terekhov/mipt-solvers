#include "csrmatrix.h"

int main(int argc, char ** argv)
{
	if( argc < 3 )
			std::cout << argv[0] << " A.mtx B.mtx [tol=1.0e-5]" << std::endl;
	else
	{
		double tol = 1.0e-5;
		std::vector<double> A, B;
		LoadVector(std::string(argv[1]),A);
		LoadVector(std::string(argv[2]),B);
		if( argc > 3 )  tol = atof(argv[3]);
		if( A.size() == B.size() ) 
		{
			double r = 0.0;
#pragma omp parallel for reduction(+:r)
			for(int i = 0; i < A.size(); ++i)
			{
				double c = A[i] - B[i];
				r += c*c;
			}
			r = std::sqrt(r);
			std::cout << "difference norm: " << r << ", a norm " << Norm(A) << ", b norm " << Norm(B) << "." << std::endl;
			if( r < tol )
			{
				std::cout << "vectors are similar" << std::endl;
				return 0;
			}
		}
		else std::cout << "vectors are not similar" << std::endl;
	}
	return 1;

}
