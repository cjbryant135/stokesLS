#ifndef __UM2HEAP_H__
#define __UM2HEAP_H__

#include "uniformmesh2d.h"
#include "heapt.h"
#include "um2heapelt.h"

namespace levelset {
        
    class UniformMesh2D_Heap : public Heap< HeapElement2D > {
                
        UniformMesh2D* mesh;
                
    public:
                
    UniformMesh2D_Heap(UniformMesh2D* m)
        : Heap< HeapElement2D >(MinHeap), mesh(m) {;}
                
        virtual ~UniformMesh2D_Heap(void) {;}
                
        void Change(const int i, const double val, const int ikhi);
        void Insert(const int i, const int j, const double val, const int ikhi);
        //   void Insert(const HeapElement2D& he);
        void Extract(int& i, int& j, const int ikhi);
                
        void Dump(std::ostream& s = std::cout, const char val = false) const;
        void Check(const int ikhi) const;
        friend class UniformMesh2D;
                
    private:
                
        void SortFrom(unsigned int n, const int ikhi);
        inline virtual void SwapElts(const int i, const int j, const int ikheap)
        {
            Swap(mesh->midata_(mesh,data[i].Index(0), data[i].Index(1), ikheap),
                 mesh->midata_(mesh,data[j].Index(0), data[j].Index(1), ikheap));
            Swap(data[i], data[j]);
        }
        inline virtual void MoveElt(const int from, const int to, const int ikheap)
        {
            mesh->midata_(mesh,data[to].Index(0), data[to].Index(1), ikheap) = -1;
            mesh->midata_(mesh,data[from].Index(0), data[from].Index(1), ikheap) = to;
            data[to] = data[from];
        }
    };
        
}

#endif // __UM2HEAP_H__

