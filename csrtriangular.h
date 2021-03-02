#ifndef _CSRTRIANGULAR_H
#define _CSRTRIANGULAR_H


#include "csrmatrix.h"
#include <assert.h>

class CSRTriangular : public CSRMatrix
{
public:
	enum Type
	{
		UpperCSR,
		LowerCSR,
		UpperCSC,
		LowerCSC
	};
private:
	Type t;
public:
	CSRTriangular(Type t) : CSRMatrix(), t(t) {}
	CSRTriangular(Type t, const std::vector<idx_t> & ia, const std::vector<idx_t> & ja, const std::vector<double> & a) : CSRMatrix(ia,ja,a), t(t) {}
	CSRTriangular(Type t, const CSRMatrix & A) : CSRMatrix(A), t(t)  {}
	CSRTriangular(const CSRTriangular & b) : CSRMatrix(b), t(b.t) {}
	CSRTriangular & operator =(CSRTriangular const & b) {t = b.t; CSRMatrix::operator =(b); return *this;}
	bool Check() const
	{
		if( t == LowerCSR || t == UpperCSC )
		{
			for(idx_t k = 0; k < Size(); ++k)
			{
				for(idx_t j = 0; j < RowSize(k); ++j)
				{
					if( Col(k,j) > k )
						return false;
				}
			}
		}
		else if( t == UpperCSR || t == LowerCSC )
		{
			for(idx_t k = 0; k < Size(); ++k)
			{
				for(idx_t j = 0; j < RowSize(k); ++j)
				{
					if( Col(k,j) < k )
						return false;
				}
			}
		}
		return true;
	}
	//The triangular matrix is stored by rows.
	//Assumes diagonal is present in each row.
	void Solve(std::vector<double> & x) const
	{
		if( t == LowerCSR )
		{
			for (idx_t k = 0; k < Size(); k++) //iterate over L part
			{
				assert(RowSize(k) > 0);
				assert(Col(k,RowSize(k)-1) == k); //diagonal
				for(idx_t r = RowSize(k) - 1; r > 0; --r)
				{
					assert(Col(k,r-1) < x.size());
					//~ assert(k > Col(k,r-1));
					x[k] -= Val(k,r-1) * x[Col(k,r-1)];
				}
				x[k] /= Val(k,RowSize(k)-1); //L diagonal
			}
		}
		else if( t == UpperCSR )
		{
			for (idx_t k = Size(); k > 0; k--) //iterate over U part
			{
				assert(RowSize(k-1) > 0);
				assert(Col(k-1,0) == k-1); //diagonal
				for(idx_t r = 1; r < RowSize(k-1); ++r)
				{
					assert(Col(k-1,r) < x.size());
					//~ assert(k-1 < Col(k-1,r));
					x[k-1] -= Val(k-1,r) * x[Col(k-1,r)];
				}
				x[k-1] /= Val(k-1,0); //U diagonal
			}
		}
		else if( t == LowerCSC )
		{
			for (idx_t k = 0; k < Size(); ++k)
			{
				assert(RowSize(k) > 0);
				assert(Col(k,0) == k); //diagonal
				x[k] /= Val(k,0);
				for (idx_t r = 1; r < RowSize(k); ++r)
				{
					assert(Col(k,r-1) < x.size());
					//~ assert(Col(k,r) > k);
					x[Col(k,r)] -= Val(k,r) * x[k];
				}
			}
		}
		else if( t == UpperCSC )
		{
			for (idx_t k = Size(); k > 0; k--) //iterate over U part
			{
				assert(RowSize(k-1) > 0);
				assert(Col(k-1,RowSize(k-1)-1) == k); //diagonal
				x[k-1] /= Val(k-1,RowSize(k-1)-1);
				for(idx_t r = 0; r < RowSize(k-1)-1; ++r)
				{
					assert(Col(k-1,r) < x.size());
					//~ assert(Col(k-1,r) > k-1);
					x[Col(k-1,r)] -= Val(k-1,r) * x[k-1];
				}
			}
		}
	}
	CSRTriangular Transpose(bool square = true) const
	{
		Type tnew;
		if( t == UpperCSR ) tnew = LowerCSR;
		else if( t == LowerCSR ) tnew = UpperCSR;
		else if( t == UpperCSC ) tnew = LowerCSC;
		else if( t == LowerCSC ) tnew = UpperCSC;
		return CSRTriangular(tnew, CSRMatrix::Transpose(square));
	}
};

#endif //_CSRTRIANGULAR_H
