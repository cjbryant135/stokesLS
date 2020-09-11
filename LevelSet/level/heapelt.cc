#include "heapelt.h"

namespace levelset {
	
template <int N>
HeapElement<N>::HeapElement(const int i, const int j, const double val) 
      : value(val)
{ ind[0] = i; ind[1] = j;}

template <int N>
HeapElement<N>::HeapElement(const int i, const int j, const int k,
                                const double val) : value(val)
{ ind[0] = i; ind[1] = j; ind[2] = k;}

template <int N>
inline int operator<(const HeapElement<N>& a, const HeapElement<N>& b)
{ return a.value < b.value;}

template <int N>
inline int operator>(const HeapElement<N>& a, const HeapElement<N>& b)
{ return a.value > b.value;}

template <int N>
inline int operator<=(const HeapElement<N>& a, const HeapElement<N>& b)
{ return a.value <= b.value;}

template <int N>
inline int operator>=(const HeapElement<N>& a, const HeapElement<N>& b)
{ return a.value >= b.value;}

}






