//
//  PM2_Boundary.h
//  LevelSet
//
//  Created by David Chopp on 7/23/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#ifndef PM2_Boundary_hpp
#define PM2_Boundary_hpp

#include "boundary2d.h"
#include "polarmesh2d.h"

namespace levelset {
	 #ifndef __BCDIR_DEFINED__
	 #define __BCDIR_DEFINED__
    enum BCDir {Right=0, Left, Up, Down};
    #endif

    class PM2_Boundary : public Boundary2D
    {
    protected:
        
        PolarMesh2D* mesh;
        
    public:
        
        PM2_Boundary(const int wxlo = 1, const int wxhi = 1,
                     const int wylo = 1, const int wyhi = 1)
        : Boundary2D(wxlo, wxhi, wylo, wyhi) {}
        
        void SetMesh(PolarMesh2D* m) {mesh = m;}
        
        void SetIBoundary(const int k, const int n);
        
        int Nabor(const int i, const int j, const BCDir d, int& ni, int& nj);
    };
    
}
#endif /* PM2_Boundary_hpp */
