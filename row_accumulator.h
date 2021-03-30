
#ifndef _ROW_ACCUMULATOR_H
#define _ROW_ACCUMULATOR_H
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include "datatype.h"

template<typename KeyType>
class RowAccumulator
{
	const idx_t NOE;
	const idx_t EOL;
	idx_t first, nnz;
	std::vector<idx_t> next;
	std::vector<KeyType> value;
public:
	RowAccumulator(idx_t num_elements, const KeyType & zero = KeyType())
	:NOE(std::numeric_limits<idx_t>::max()), EOL(NOE-1)
	{
		first = EOL;
		nnz = 0;
		next.resize(num_elements,NOE);
		value.resize(num_elements,zero);
	}
	idx_t Begin() const {return first;}
	idx_t End() const {return EOL;}
	idx_t Next(idx_t i) const {return next[i];}
	KeyType & operator [](idx_t i)
	{
		if( next[i] == NOE )
		{
			next[i] = first;
			first = i;
			nnz++;
		}
		return value[i];
	}
	KeyType & Get(idx_t i) {return value[i];}
	const KeyType & Get(idx_t i) const {return value[i];}
	bool Contains(idx_t i) const {return next[i] != NOE;}
	void Put(idx_t i1, idx_t i2, const std::vector<idx_t> & ja, const std::vector<KeyType> & a, const KeyType & coef)
	{
		assert(Empty());
		for(idx_t i = i1; i < i2; ++i)
		{
			value[ja[i]] = a[i]*coef;
			next[ja[i]] = first;
			first = ja[i];
			nnz++;
		}
	}
	void Add(idx_t i1, idx_t i2, const std::vector<idx_t> & ja, const std::vector<KeyType> & a, const KeyType & coef)
	{
		for(idx_t i = i1; i < i2; ++i)
		{
			value[ja[i]] += a[i]*coef;
			if( next[ja[i]] == NOE )
			{
				next[ja[i]] = first;
				first = ja[i];
				nnz++;
			}
		}
	}
	//assumes list is sorted and insertions are sorted
	idx_t InsertOrdered(idx_t curr, idx_t ind, const KeyType & val)
	{
		if (Contains(ind))
		{
			Get(ind) += val;
			return ind;
		}
		else if (1 + val != 1)
		{
			if (first == EOL || first > ind) //no entries
			{
				next[ind] = first;
				first = ind;
				value[ind] = val;
			}
			else
			{
				idx_t foll = next[curr] == NOE ? first : curr;
				while (foll < ind)
				{
					curr = foll;
					foll = next[curr];
				}
				next[curr] = ind;
				next[ind] = foll;
				value[ind] = val;
			}
			return ind;
		}
		else return curr;
	}
	void Sort()
	{
		std::vector<idx_t> cols(nnz);
		idx_t i = first, k = 0;
		while(i != EOL)
		{
			cols[k++] = i;
			i = next[i];
		}
		if( !cols.empty() )
		{
			std::sort(cols.begin(),cols.end());
			first = cols[0];
			for(idx_t q = 1; q < (idx_t)cols.size(); ++q)
				next[cols[q-1]] = cols[q];
			next[cols.back()] = EOL;
		}
	}
	void Mult(const KeyType & coef)
	{
		idx_t i = first;
		while( i != EOL )
		{
			value[i] *= coef;
			i = next[i];
		}
	}
	void Get(std::vector<idx_t> & ja, std::vector<KeyType> & a) const
	{
		idx_t i = first;
		while( i != EOL )
		{
			if( 1+value[i] != 1 )
			{
				ja.push_back(i);
				a.push_back(value[i]);
			}
			i = next[i];
		}
	}
	void Clear(const KeyType & zero = KeyType())
	{
		idx_t i = first, j;
		while(i != EOL)
		{
			j = next[i];
			value[i] = zero;
			next[i] = NOE;
			i = j;
		}
		first = EOL;
		nnz = 0;
	}
	void Print() const
	{
		idx_t i = first;
		while( i != EOL )
		{
			std::cout << "[" << i << "]: " << value[i] << std::endl;
			i = next[i];
		}
	}
	idx_t Size() const {return nnz;}
	bool Empty() const {return first == EOL;}
};

#endif //_LINKED_LIST_H
