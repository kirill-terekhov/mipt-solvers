#ifndef _CONDEST_H
#define _CONDEST_H

#include "csrmatrix.h"
#include "row_accumulator.h"
#include <cmath>

/*
 * Invnorm - Running estimator of inverse norms for triangular factors in LU methods
 * Implemented following
 * [1] "A robust and efficient ILU that incorporates the growth of the inverse triangular factors"
 * by M.Bollhoefer
 * [2] "A robust ILU with pivoting based on monitoring the growth of the inverse factors"
  * by M.Bollhoefer
 */

class Invnorm
{
	std::vector<double> eta, zeta;
	double c_eta, c_zeta;
public:
	Invnorm(idx_t Size) : eta(Size,0.0), zeta(Size,0.0) {}
	Invnorm(const Invnorm & b) :eta(b.eta), zeta(b.zeta) {}
	Invnorm & operator = (Invnorm const & b) {eta = b.eta; zeta = b.zeta; return *this;}
	double Estimate(const RowAccumulator<double> & row, idx_t k)
	{
		double mup, mun, nup, nun, sp, sn, c1, c2;
		mup = -eta[k] + 1.0;
		mun = -eta[k] - 1.0;
		sp = sn = 0;
		for(idx_t j = row.Begin(); j != row.End(); j = row.Next(j)) if( j > k )
		{
			sp += std::fabs(eta[j]+row.Get(j)*mup);
			sn += std::fabs(eta[j]+row.Get(j)*mun); 
		}
		if( std::fabs(mup) + sp > std::fabs(mun) + sn )
			c_eta = mup;
		else 
			c_eta = mun;
		c1 = std::max(std::fabs(mup),std::fabs(mun));
		mup = -zeta[k] + 1.0;
		mun = -zeta[k] - 1.0;
		sp = sn = 0;
		for(idx_t j = row.Begin(); j != row.End(); j = row.Next(j)) if( j > k )
		{
			nup = std::fabs(zeta[j]+row.Get(j)*mup);
			nun = std::fabs(zeta[j]+row.Get(j)*mun);
			if( nup > std::max(2.0*std::fabs(zeta[j]),0.5) ) sp++;
			if( nun > std::max(2.0*std::fabs(zeta[j]),0.5) ) sn++;
			if( std::fabs(zeta[j]) > std::max(2*nup,0.5) ) sp--;
			if( std::fabs(zeta[j]) > std::max(2*nun,0.5) ) sn--;
		}
		if( std::fabs(mup) + sp > std::fabs(mun) + sn )
			c_zeta = mup;
		else
			c_zeta = mun;
		c2 = std::max(std::fabs(mup),std::fabs(mun));
		return std::max(c1,c2);
	}
	void Update(const RowAccumulator<double> & row, idx_t k)
	{
		for(idx_t j = row.Begin(); j != row.End(); j = row.Next(j)) if( j > k )
		{
			eta[j]  += row.Get(j)*c_eta;
			zeta[j] += row.Get(j)*c_zeta;
		}
	}
	size_t Bytes() {return get_bytes(eta) + get_bytes(zeta);}
};

#endif //_CONDEST_H
