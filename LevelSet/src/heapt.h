#ifndef __HEAP_H__
#define __HEAP_H__

#include "datablock.h"

namespace levelset {
        
    template <class T>
        inline void Swap(T& a, T& b) {T t = a; a = b; b = t;}
        
    enum HeapType {MinHeap=0, MaxHeap=1};
        
    template <class T>
        class Heap {
    protected:
                
        DataBlock<T>  data;
        HeapType      type;
        unsigned int  length;
                
    public:
                
    Heap(HeapType h=MinHeap) : data(), type(h), length(0) {;}
                
        void Insert(const T& elt);
        T Extract(void);
        inline T Top(void) const {return data[0];}
                
        inline char Empty(void) {return length==0;}
                
#ifdef LEVEL_DEBUG
        void DumpHeap(const int n) const;
        void Check(void) const;
        friend class UniformMesh2D;
#endif
                
    protected:
                
        enum HeapChild {HLeft=0, HRight=1};
                
        void SortFrom(unsigned int n);
        inline void SwapElts(const int i, const int j) 
        { Swap(data[i], data[j]); }
        inline void MoveElt(const int from, const int to)
        {data[to] = data[from];} 
                
        inline unsigned int parent(const int k) const {return (k-1)/2;}
        inline unsigned int child(const int k, const HeapChild i) const 
        {return 2*k+i+1;}
                
    };
        
    inline int heapT_log2(const unsigned int n)
    {
        int l,k; 
        for (l=0, k=n; k>0; ++l, k/=2);
        return l;
    }
        
    inline int heapT_row2(const unsigned int n)
    {
        return (0x01 << heapT_log2(n)) - 1;
    }
        
}

#ifndef NO_INCLUDE_CC_FILE
#include "heapt.cc"
#endif // NO_INCLUDE_CC_FILE

#endif // __HEAP_H__

