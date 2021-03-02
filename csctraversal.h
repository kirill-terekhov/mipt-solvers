#ifndef _CSCTRAVERSAL_H
#define _CSCTRAVERSAL_H
#include "csrmatrix.h"
#include <assert.h>
/*
 * Traverse transposed matrix without assembly
 * only square matrices
 */
class CSCTraversal
{
	const idx_t NOE;
	const idx_t EOL;
	std::vector<idx_t> curr;
	std::vector<idx_t> first;
	std::vector<idx_t> next;
	const CSRMatrix & A;
	idx_t col;
public:
	CSCTraversal(const CSRMatrix & A, idx_t Size) 
	: NOE(std::numeric_limits<idx_t>::max()), EOL(NOE-1), 
	  curr(Size,NOE), first(Size,EOL), next(Size,NOE), A(A)
	{
		for(idx_t k = A.Size(); k > 0; --k)
			NewRow(k-1);
		col = 0;
	}
	idx_t Begin() const {return first[col];}
	idx_t Next(idx_t i) const {return next[i];}
	idx_t End() const {return EOL;}
	idx_t Position(idx_t i) const {return curr[i];}
	idx_t CurrentColumn() const {return col;}
	void NewRow(idx_t i)
	{
		if ( A.RowSize(i) > 0)
		{
			curr[i] = 0;
			idx_t j = A.Col(i,curr[i]);
			assert(i != first[j]);
			next[i] = first[j];
			first[j] = i;
		}
	}
	void NextColumn()
	{
		idx_t k, l, j, cur, nxt;
		k = first[col++];
		while(k != EOL)
		{
			l = next[k];
			curr[k]++;
			if( A.RowSize(k) - curr[k] > 0 )
			{
				j = A.Col(k,curr[k]);
				//ordered insert to linked-list
				if( first[j] > k )
				{
					next[k] = first[j];
					first[j] = k;
				}
				else
				{
					cur = nxt = first[j];
					while( nxt < k )
					{
						cur = nxt;
						nxt = next[cur];
					}
					assert(cur < k && k < nxt);
					next[cur] = k;
					next[k] = nxt;
				}
			}
			else curr[k] = NOE;
			k = l;
		}
	}
	size_t Bytes() const {return get_bytes(curr) + get_bytes(first) + get_bytes(next) + 2*sizeof(int) + sizeof(const CSRMatrix &);}
};

#endif //_CSCTRAVERSAL_H
