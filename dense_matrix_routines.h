#ifndef _DENSE_MATRIX_ROUTINES_H
#define _DENSE_MATRIX_ROUTINES_H
#include "priority_queue.h"
//   n        
// -----
// | A | m
// -----
static void printmxn(const double* a, int m, int n)
{
	std::cout << "matrix(";
	for (int l = 0; l < m; ++l)
	{
		std::cout << std::endl << "\t[";
		for (int q = 0; q < n; ++q)
			std::cout << std::setw(12) << a[l * n + q] << ",";
		std::cout << "\b],";
	}
	std::cout << "\b)" << std::endl;
}
//   n          m
// -----      -----
// | A | m -> | B | n
// -----      -----
static void transposemxn(const double* a, double* b, int m, int n, const int* ord = NULL)
{
	if (ord)
	{
		for (int j = 0; j < m; j++)
			for (int i = 0; i < n; i++)
				b[i * m + j] = a[j * n + ord[i]];
	}
	else
	{
		for (int j = 0; j < m; j++)
			for (int i = 0; i < n; i++)
				b[i * m + j] = a[j * n + i];
	}
}
//   n          n
// -----      -----
// | A | n -> | B | n
// -----      -----
static void transposenxn(double* a, int n)
{
	for (int j = 0; j < n; j++)
		for (int i = j + 1; i < n; i++)
			std::swap(a[i * n + j], a[j * n + i]);
}
//   k         n          n
// -----     -----      -----
// | A | m * | B | k =  | C | m
// -----     -----      -----
static void matmulmxnxk(const double* a, const double* b, double* c, int m, int n, int k)
{
	int i, j, l;
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			c[j * n + i] = 0.0;
			for (l = 0; l < k; l++)
				c[j * n + i] += a[j * k + l] * b[l * n + i];
		}
}
// note: A is row major, but b and x are column-major
static int solvenxnxm(double * A, double * x, double * b, int n, int m, int * order)
{
	double max;
	for(int i = 0; i < n; i++) 
		order[i] = i;
		
	for(int i = 0; i < n; i++)
	{
		int maxk = i, maxq = i;
		max = fabs(A[maxk*n+maxq]);
		//Find best pivot
		for(int q = i; q < n; q++) // over columns
			for(int k = i; k < n; k++) // over rows
				if( fabs(A[k*n+q]) > max )
				{
					max = fabs(A[k*n+q]);
					maxk = k;
					maxq = q;
				}
		//Exchange rows
		if( maxk != i )
		{
			for(int q = 0; q < n; q++)
				std::swap(A[maxk * n + q], A[i * n + q]);
			//exchange rhs
			for(int k = 0; k < m; k++)
				std::swap(b[maxk + k * n], b[i + k * n]);
		}
		//Exchange columns
		if( maxq != i )
		{
			for(int k = 0; k < n; k++)
				std::swap(A[k * n + maxq], A[k * n + i]);
			//remember order in sol
			std::swap(order[maxq], order[i]);
		}
			
		if( 1 + A[i*n+i] == 1 )
			return i+1;
		for(int k = i+1; k < n; k++)
		{
			A[i*n+k] /= A[i*n+i];
			A[k*n+i] /= A[i*n+i];
		}
		for(int k = i+1; k < n; k++)
			for(int q = i+1; q < n; q++)
				A[k*n+q] -= A[k*n+i] * A[i*n+i] * A[i*n+q];
		for(int k = 0; k < m; k++)
		{
			for(int j = i+1; j < n; j++) //iterate over columns of L
				b[j+k*n] -= b[i+k*n] * A[j*n+i];
			b[i+k*n] /= A[i*n+i];
		}
	}
		
	for(int k = 0; k < m; k++)
	{
		for(int i = n-1; i >= 0; i--) //iterate over rows of U
			for(int j = i+1; j < n; j++)
				b[i+k*n] -= b[j+k*n] * A[i*n+j];
		for(int i = 0; i < n; i++)
			x[order[i]+k*n] = b[i+k*n];
	}
	return 0;
}


static void unitnxn(double* u, int n)
{
	std::fill(u, u + n * n, 0.0);
	for (int l = 0; l < n; ++l) u[l * n + l] = 1.0;
}

#endif //_DENSE_MATRIX_ROUTINES_H
