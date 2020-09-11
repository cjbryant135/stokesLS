//
//  pm2heapelt.h
//  LevelSet
//
//  Created by David Chopp on 7/25/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#ifndef pm2heapelt_h
#define pm2heapelt_h

#include "heapelt2d.h"
#include "polarmesh2d.h"

namespace levelset {
    
    class PM2_HeapElement : public HeapElement2D {
        
        PolarMesh2D* mesh;
        int ikheap;
        
    public:
        
        PM2_HeapElement(PolarMesh2D* m, const int i, const int j,
                        const int iktemp, const double val)
        : HeapElement2D(i,j,val), mesh(m), ikheap(iktemp) {;}
        PM2_HeapElement(void) : HeapElement2D(), mesh(NULL) {;}
        
        inline PM2_HeapElement& operator=(const PM2_HeapElement& a)
        {mesh->midata_(mesh,ind[0],ind[1],ikheap)
            = a.mesh->midata_(a.mesh,a.ind[0],a.ind[1],ikheap);
            HeapElement2D::operator=(a);
            return *this;}
        
#ifdef LEVEL_DEBUG
        friend class PolarMesh2D;
#endif
    };
    
}
#endif /* pm2heapelt_h */
