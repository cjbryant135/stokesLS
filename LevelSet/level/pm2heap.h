//
//  pm2heap.hpp
//  LevelSet
//
//  Created by David Chopp on 7/25/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#ifndef pm2heap_hpp
#define pm2heap_hpp

#include "polarmesh2d.h"
#include "heapt.h"
#include "pm2heapelt.h"

namespace levelset {
    
    class PolarMesh2D_Heap : public Heap< HeapElement2D > {
        
        PolarMesh2D* mesh;
        
    public:
        
        PolarMesh2D_Heap(PolarMesh2D* m)
        : Heap< HeapElement2D >(MinHeap), mesh(m) {;}
        
        virtual ~PolarMesh2D_Heap(void) {;}
        
        void Change(const int i, const double val, const int ikhi);
        void Insert(const int i, const int j, const double val, const int ikhi);
        //   void Insert(const HeapElement2D& he);
        void Extract(int& i, int& j, const int ikhi);
        
        void Dump(std::ostream& s = std::cout, const char val = false) const;
        void Check(const int ikhi) const;
        friend class PolarMesh2D;
        
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
#endif /* pm2heap_hpp */
