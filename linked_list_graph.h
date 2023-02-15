#ifndef _LINKED_LIST_GRAPH_H
#define _LINKED_LIST_GRAPH_H

#include "datatype.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <assert.h>

class LinkedListGraph
{
	const idx_t NOE;
	const idx_t EOL;
	idx_t first, nnz;
	std::vector<idx_t> next;
public:
	LinkedListGraph(idx_t num_elements)
	:NOE(std::numeric_limits<idx_t>::max()), EOL(NOE-1)
	{
		first = EOL;
		nnz = 0;
		next.resize(num_elements,NOE);
	}
	idx_t Begin() const {return first;}
	idx_t End() const {return EOL;}
	idx_t Next(idx_t i) const {return next[i];}
	bool Contains(idx_t i) const {return next[i] != NOE;}
	void Put(idx_t i1, idx_t i2, const std::vector<idx_t> & ja)
	{
		assert(Empty());
		for(idx_t i = i1; i < i2; ++i)
		{
			next[ja[i]] = first;
			first = ja[i];
			nnz++;
		}
	}
	void Add(idx_t i1, idx_t i2, const std::vector<idx_t> & ja)
	{
		for(idx_t i = i1; i < i2; ++i)
		{
			if( next[ja[i]] == NOE )
			{
				next[ja[i]] = first;
				first = ja[i];
				nnz++;
			}
		}
	}
	void Set(idx_t i)
	{
		if (next[i] == NOE)
		{
			next[i] = first;
			first = i;
			nnz++;
		}
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
			for(idx_t q = 1; q < cols.size(); ++q)
				next[cols[q-1]] = cols[q];
			next[cols.back()] = EOL;
		}
	}
	void Get(std::vector<idx_t> & ja) const
	{
		idx_t i = first;
		while( i != EOL )
		{
			ja.push_back(i);
			i = next[i];
		}
	}
	void Clear()
	{
		idx_t i = first, j;
		while(i != EOL)
		{
			j = next[i];
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
			std::cout << i << " ";
			i = next[i];
		}
		std::cout << std::endl;
	}
	idx_t Size() const {return nnz;}
	bool Empty() const {return first == EOL;}
};

#endif //_LINKED_LIST_GRAPH_H
