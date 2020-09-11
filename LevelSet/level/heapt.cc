#include <iostream>
#include <iomanip>
#include <stdlib.h>
//#include "heapt.h"

namespace levelset {
	
	template <class T>
	void Heap<T>::Insert(const T& elt)
	{
		++length;
		if (length > data.Length()) {
			data.Resize(heapT_row2(length));
		}
		data[length-1] = elt;
		SortFrom(length-1);
	}
	
	template <class T>
	T Heap<T>::Extract(void)
	{
		T top = data[0];
		--length;
		if (length > 0) {
			MoveElt(length,0);
			SortFrom(0);
		}
		return top;
	}
	
	
	template <class T>
	void Heap<T>::SortFrom(unsigned int n)
	{
		if (type == MinHeap) {
			int k;
			for (k=n; k>0; k=parent(k)) 
				if (data[k] < data[parent(k)]) 
					SwapElts(k,parent(k));
				else 
					break;
			
			int c;
			while (child(k,HLeft) < length) {
				if (child(k,HRight) >= length) 
					c = child(k,HLeft);
				else
					c = (data[child(k,HLeft)] < data[child(k,HRight)]) ? 
					child(k,HLeft) : child(k,HRight);
				if (data[k] > data[c]) {
					SwapElts(k,c);
					k = c;
				} else 
					break;
			}
		} else {
			int k;
			for (k=n; k>0; k=parent(k)) 
				if (data[k] > data[parent(k)]) 
					SwapElts(k,parent(k));
				else 
					break;
			
			int c;
			while (child(k,HLeft) < length) {
				if (child(k,HRight) >= length) 
					c = child(k,HLeft);
				else
					c = (data[child(k,HLeft)] > data[child(k,HRight)]) ? 
					child(k,HLeft) : child(k,HRight);
				if (data[k] < data[c]) {
					SwapElts(k,c);
					k = c;
				} else 
					break;
			}
		}
		
	}
	
#ifdef LEVEL_DEBUG
	template <class T>
	void Heap<T>::DumpHeap(const int w) const
	{
		int depth = heapT_log2(length);
		int rows = 0x01 << (depth-1);
		
		std::cout << std::setfill(' ');
		for (int row=0; row<rows; ++row) {
			for (int d=0; d<depth; ++d) 
				if (!(row%(0x01 << (depth-d-1)))) {
					int b = (0x01 << d) - 1;
					int r = row/(0x01 << (depth-d-1));
					if (b+r < length)
						std::cout << std::setw(w) << data[b+r];
					else
						std::cout << std::setw(w) << " ";
				} else 
					std::cout << std::setw(w) << " ";
				std::cout << "\n";
		}
	}
	
	template <class T>
	void Heap<T>::Check(void) const
	{
		if (type == MinHeap) {
			for (int ii=0; ii<length-1; ++ii) {
				if (child(ii,HLeft) < length) {
					if (data[ii] > data[child(ii,HLeft)]) {
						std::cerr << "Error in Heap at element " << ii << " to the left\n";
						exit(1);
					}
				}
				if (child(ii,HRight) < length) {
					if (data[ii] > data[child(ii,HRight)]) {
						std::cerr << "Error in Heap at element " << ii << " to the right\n";
						exit(1);
					}
				}
			}
		}
		else {
			for (int ii=0; ii<length-1; ++ii) {
				if (child(ii,HLeft) < length) {
					if (data[ii] < data[child(ii,HLeft)]) {
						std::cerr << "Error in Heap at element " << ii << " to the left\n";
						exit(1);
					}
				}
				if (child(ii,HRight) < length) {
					if (data[ii] < data[child(ii,HRight)]) {
						std::cerr << "Error in Heap at element " << ii << " to the right\n";
						exit(1);
					}
				}
			}
		}
	}
	
#endif
}

