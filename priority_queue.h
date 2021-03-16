#ifndef _PRIORITY_QUEUE_H
#define _PRIORITY_QUEUE_H
#include <vector>
#include <iostream>
#include "datatype.h"

template<typename KeyType, typename Comparator = std::less<KeyType> >
class PriorityQueue
{
	const idx_t NOE;
	idx_t size_max, size_cur;
	std::vector<idx_t> heap;
	std::vector<idx_t> index;
	std::vector<KeyType> keys;
	void Swap(idx_t i, idx_t j)
	{
		idx_t t = heap[i];
		heap[i] = heap[j];
		heap[j] = t;
		index[heap[i]] = i;
		index[heap[j]] = j;
	}
	void BubbleUp(idx_t k)
	{
		while(k > 1 && !Comparator()(keys[heap[k/2]],keys[heap[k]]))
		{
			Swap(k, k/2);
			k = k/2;
		}
	}
	void BubbleDown(idx_t k)
	{
		idx_t j;
		while(2*k <= size_cur)
		{
			j = 2*k;
			if(j < size_cur && !Comparator()(keys[heap[j]],keys[heap[j+1]]) )
				j++;
			if(Comparator()(keys[heap[k]],keys[heap[j]]) )
				break;
			Swap(k, j);
			k = j;
		}
	}
public:
	PriorityQueue(idx_t num_elements, const KeyType & c = KeyType())
	:NOE(std::numeric_limits<idx_t>::max())
	{
		size_max = num_elements;
		size_cur = 0;
		keys.resize(size_max,c);
		heap.resize(size_max+1,NOE);
		index.resize(size_max+1,NOE);
	}
	void Push(idx_t i, const KeyType & key)
	{
		size_cur++;
		index[i] = size_cur;
		heap[size_cur] = i;
		keys[i] = key;
		BubbleUp(size_cur);
	}
	idx_t Pop()
	{
		if( size_cur == 0 ) 
			return NOE;
		idx_t min = heap[1];
		Swap(1, size_cur--);
		BubbleDown(1);
		index[min] = NOE;
		heap[size_cur+1] = NOE;
		return min;
	}
	idx_t Peek() const
	{
		if( size_cur == 0 ) return NOE;
		return heap[1];
	}
	void DecreaseKey(idx_t i, const KeyType & key)
	{
		keys[i] = key;
		BubbleUp(index[i]);
	}
	void IncreaseKey(idx_t i, const KeyType & key)
	{
		keys[i] = key;
		BubbleDown(index[i]);
	}
	void ChangeKey(idx_t i, const KeyType & key)
	{
		keys[i] = key;
		if( !Comparator()(keys[i],key) )	
			BubbleUp(index[i]);
		else
			BubbleDown(index[i]);
	}
	void ResetKey(idx_t i, const KeyType & key) { keys[i] = key; }
	const KeyType & GetKey(idx_t i) const {return keys[i];}
	void Clear() {while( !Empty() ) Pop();}
	bool Contains(idx_t i) const {return index[i] != NOE;}
	idx_t  Size() const {return size_cur;}
	bool Empty() const {return size_cur == 0;}
};

#endif //_PRIORITY_QUEUE_H
