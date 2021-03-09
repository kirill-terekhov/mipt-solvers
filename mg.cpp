#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstring>
#include <vector>
#include <iomanip>
//adopted from https://people.sc.fsu.edu/~jburkardt/cpp_src/multigrid_poisson_1d/multigrid_poisson_1d.cpp

// transfers data from a coarse to a finer grid.
void ctof ( const std::vector<double> & uc, // coarse error
			std::vector<double> & uf ) // fine unknowns
{
	for (size_t ic = 0; ic < uc.size(); ic++)
		uf[ic*2] += uc[ic];
	for (size_t ic = 0; ic < uc.size() - 1; ic++)
		uf[ic*2+1] += 0.5 * ( uc[ic] + uc[ic+1] );
}
// transfers data from a fine grid to a coarser grid.
void ftoc ( const std::vector<double> & uf, // fine unknowns
			const std::vector<double> & rf, // fine rhs
			std::vector<double> & uc, // coarse correction (init to zero)
			std::vector<double> & rc) // coarse rhs
{
	std::fill(uc.begin(),uc.end(),0);
	rc[0] = 0.0;
	for (size_t ic = 1; ic < rc.size() - 1; ic++)
		rc[ic] = 4.0 * ( rf[ic*2] + uf[ic*2-1] - 2.0 * uf[ic*2] + uf[ic*2+1] );
	rc[rc.size()-1] = 0.0;
}

// carries out one step of a Gauss-Seidel iteration.
double gauss_seidel(const std::vector<double> & r, //rhs
					std::vector<double> & u ) //solution estimation
{
	double norm = 0.0;
	for (size_t i = 1; i < u.size() - 1; i++ )
	{
		double u_old = u[i];
		u[i] = 0.5 * ( u[i-1] + u[i+1] + r[i] );
		norm += std::fabs( u[i] - u_old );
	}
	return norm;
}

// solve a 1d poission problem using gauss-seidel method
void singlegrid_poisson_1d ( 
	size_t n, //number of cells
	double a, //left-most node position
	double b, //right-most node position
	double ua, //left boundary condition
	double ub, //right boundary condition
	double force ( double x ), //external force
	double exact ( double x )) //reference solution
{
	const double tol = 1.0e-3;
	std::vector<double> r(n+1); //right hand side
	std::vector<double> u(n+1); //solution
	double h = (b - a)/static_cast<double>(n);
	r[0] = ua;
	for(size_t i = 0; i < n+1; ++i)
	{
		double x = a + h*i;
		r[i] = h*h + force(x);
	}
	r[n] = ub;
	std::fill(u.begin(),u.end(),0.0);
	
	size_t it_num = 0;
	double diff = 1.0e+20;
	while( diff > tol ) 
	{
		diff = gauss_seidel(r,u);
		std::cout << "iteration " << std::setw(10) << it_num << " error " << std::setw(14) << diff << "\r";
		std::cout.flush();
		it_num++;
	}
	std::cout << std::endl;
	
	double err = 0;
	for(size_t i = 0; i < n+1; ++i)
	{
		double x = a + h*i;
		err += std::pow(exact(x) - u[i]*h*h,2);
		//std::cout << "at " << std::setw(14) << x << " exact " << exact(x) << " obtained " << u[i]*h*h << std::endl;
	}
	err = std::sqrt(err);
	
	std::cout << "Converged in " << it_num << " Gauss-Seidel iterations below " << tol << " tolerance, true solution error " << err << std::endl;	
}

//solve a 1d poisson problem using multigrid method
void multigrid_poisson_1d ( 
	size_t n, //number of cells
	double a, //left-most node position
	double b, //right-most node position
	double ua, //left boundary condition
	double ub, //right boundary condition
	double force ( double x ), //external force
	double exact ( double x )) //reference solution
{
	const double tol = 1.0e-3, utol = 0.7;
	size_t k = std::log2 ( n );
	if ( n != std::pow ( 2, k ) ) 
	{
		std::cout << "Number of cells should be power of 2, got " << n << std::endl;
		return;
	}
	std::vector< std::vector<double> > r(k); //right hand side for k levels
	std::vector< std::vector<double> > u(k); //solution for k levels
	for(size_t i = 0; i < k; ++i) //initialize each level
	{
		r[i].resize( n / std::pow( 2, i ) + 1, 0.0);
		u[i].resize( n / std::pow( 2, i ) + 1, 0.0);
	}
	//initialize first level rhs
	double h = (b - a)/static_cast<double>(n);
	r[0][0] = ua;
	for(size_t i = 0; i < n+1; ++i)
	{
		double x = a + h*i;
		r[0][i] = h*h + force(x);
	}
	r[0][n] = ub;
	size_t it_num = 0;
	size_t level = 0;
	double diff0, diff1 = 0;
	int j = 0;
	while( true )
	{
		diff0 = diff1;
		diff1 = gauss_seidel(r[level],u[level]);
		std::cout << "iteration " << std::setw(10) << it_num << " level " << std::setw(2) << level << " error " << std::setw(14) << diff1 << "\r";
		std::cout.flush();
		it_num++;
		j++;
		if( j < 4 ) //at least 4 gauss-seidel iterations
			continue;
		else if( diff1 < tol && level == 0 ) //satisfactory convergence on finest level - exit
			break;
		else if( diff1 < tol ) //satisfactory convergence on coarse level, go finer
		{
			ctof(u[level],u[level-1]);
			level--;
			j = 0;
		}
		else if( utol * diff0 <= diff1 && level < k-1 ) // enough iterations, go coarser
		{
			ftoc(u[level],r[level],u[level+1],r[level+1]);
			level++;
			j = 0;
		}
	}
	std::cout << std::endl;
	
	double err = 0;
	for(size_t i = 0; i < n+1; ++i)
	{
		double x = a + h*i;
		err += std::pow(exact(x) - u[0][i]*h*h,2);
		//std::cout << "at " << std::setw(14) << x << " exact " << exact(x) << " obtained " << u[0][i]*h*h << std::endl;
	}
	err = std::sqrt(err);
	
	
	std::cout << "Converged in " << it_num << " Multigrid Gauss-Seidel iterations below " << tol << " tolerance, true solution error " << err << std::endl;	
	
	return;
}

double force( double x )
{
  return - x * ( x + 3.0 ) * std::exp ( x );
}

double exact( double x )
{
  return x * ( x - 1.0 ) * std::exp ( x );
}

int main()
{
	const int N = 512;
	singlegrid_poisson_1d(N,0,1,0.0,0.0,force,exact);
	multigrid_poisson_1d(N,0,1,0.0,0.0,force,exact);
	return 0;
}



