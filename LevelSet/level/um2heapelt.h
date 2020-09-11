#ifndef __UM2HEAPELT_H__
#define __UM2HEAPELT_H__

#include "heapelt2d.h"
#include "uniformmesh2d.h"

namespace levelset {
        
    class UM2_HeapElement : public HeapElement2D {

        UniformMesh2D* mesh;
        int ikheap;
   
    public:

    UM2_HeapElement(UniformMesh2D* m, const int i, const int j,
                    const int iktemp, const double val)
        : HeapElement2D(i,j,val), mesh(m), ikheap(iktemp) {;}
    UM2_HeapElement(void) : HeapElement2D(), mesh(NULL) {;}
   
        inline UM2_HeapElement& operator=(const UM2_HeapElement& a)
            {mesh->midata_(mesh,ind[0],ind[1],ikheap)
             = a.mesh->midata_(a.mesh,a.ind[0],a.ind[1],ikheap);
             HeapElement2D::operator=(a);
             return *this;}

#ifdef LEVEL_DEBUG
        friend class UniformMesh2D;
#endif
    };

}               
#endif // __HEAPELT_H__
