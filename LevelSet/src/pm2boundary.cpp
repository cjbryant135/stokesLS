//
//  PM2_Boundary.cpp
//  LevelSet
//
//  Created by David Chopp on 7/23/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#include "pm2boundary.h"

namespace levelset {
    
    void PM2_Boundary::SetIBoundary(const int k, const int n)
    {
        int i, j;
        
        for (i=-width[BDRY2D_XLO]; i<0; ++i)
            for (j=0; j<mesh->maxj; ++j)
                mesh->midata_(mesh,i,j,k) = n;
        
        for (i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i)
            for (j=0; j<mesh->maxj; ++j)
                mesh->midata_(mesh,i,j,k) = n;
        
        for (j=-width[BDRY2D_YLO]; j<0; ++j)
            for (i=0; i<mesh->maxi; ++i)
                mesh->midata_(mesh,i,j,k) = n;
        
        for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
            for (i=0; i<mesh->maxi; ++i)
                mesh->midata_(mesh,i,j,k) = n;
        
        for (i=-width[BDRY2D_XLO]; i<0; ++i) {
            for (j=-width[BDRY2D_YLO]; j<0; ++j)
                mesh->midata_(mesh,i,j,k) = n;
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
                mesh->midata_(mesh,i,j,k) = n;
        }
        
        for (i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i) {
            for (j=-width[BDRY2D_YLO]; j<0; ++j)
                mesh->midata_(mesh,i,j,k) = n;
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
                mesh->midata_(mesh,i,j,k) = n;
        }
    }
    
    int PM2_Boundary::Nabor(const int i, const int j, const BCDir d, int& ni, int& nj)
    {
        int valid = 1;
        switch(d) {
            case Right: if (i < mesh->maxi) {ni=i+1; nj=j;} else valid=0; break;
            case Left: if (i > 0) {ni=i-1; nj=j;} else valid=0; break;
            case Up: if (j < mesh->maxj) {ni=i; nj=j+1;} else valid=0; break;
            case Down: if (j > 0) {ni=i; nj=j-1;} else valid=0; break;
        }
        return valid;
    }
    
}
