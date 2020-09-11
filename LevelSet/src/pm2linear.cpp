//
//  pm2linear.cpp
//  LevelSet
//
//  Created by David Chopp on 7/23/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#include "pm2linear.h"

namespace levelset {
    
    void PM2_LinearBdry::Apply(const int k)
    {
        int i, j;
        
        for (i=-width[BDRY2D_XLO]; i<0; ++i)
            for (j=0; j<mesh->maxj; ++j)
                mesh->mdata_(mesh,i,j,k) = (1-i)*mesh->mdata_(mesh,0,j,k)
                +i*mesh->mdata_(mesh,1,j,k);
        
        for (i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i)
            for (j=0; j<mesh->maxj; ++j)
                mesh->mdata_(mesh,i,j,k) = (i-mesh->maxi+2)*mesh->mdata_(mesh,mesh->maxi-1,j,k)
                -(i-mesh->maxi+1)*mesh->mdata_(mesh,mesh->maxi-2,j,k);
        
        int ydiff = mesh->Index(0,mesh->maxj,0) - mesh->Index(0,0,0);
        
        for (j=0; j<width[BDRY2D_YLO]; ++j)
            for (i=0; i<mesh->tmaxi; ++i)
                mesh->thedata[mesh->Index(i,j,k)] = mesh->thedata[mesh->Index(i,j,k)+ydiff];

        for (j=0; j<width[BDRY2D_YHI]; ++j)
            for (i=0; i<mesh->tmaxi; ++i)
                mesh->thedata[mesh->Index(i,mesh->tmaxj-1-j,k)] = mesh->thedata[mesh->Index(i,mesh->tmaxj-1-j,k)-ydiff];
        
    }
    
    void PM2_LinearBdry::Apply(const int ic, const int jc, const int k)
    {
        if (ic == 0 || ic == 1)
            for (int i=-width[BDRY2D_XLO]; i < 0; ++i)
                mesh->mdata_(mesh,i, jc, k) = (1-i)*mesh->mdata_(mesh,0,jc,k)
                +i*mesh->mdata_(mesh,1,jc,k);
        if (ic == mesh->maxi-1 || ic == mesh->maxi-2)
            for (int i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i)
                mesh->mdata_(mesh,i,jc,k) = (i-mesh->maxi+2)*mesh->mdata_(mesh,mesh->maxi-1,jc,k)
                -(i-mesh->maxi+1)*mesh->mdata_(mesh,mesh->maxi-2,jc,k);
        if (ic > mesh->maxi-1-width[BDRY2D_XLO]) {
            mesh->mdata_(mesh,ic-mesh->maxi, jc, k) = mesh->mdata_(mesh,ic,jc,k);
            if (jc > mesh->maxj-1-width[BDRY2D_YLO])
                mesh->mdata_(mesh,ic-mesh->maxi, jc-mesh->maxj, k) = mesh->mdata_(mesh,ic,jc,k);
        }
        if (jc > mesh->maxj-1-width[BDRY2D_YLO]) {
            mesh->mdata_(mesh,ic, jc-mesh->maxj, k) = mesh->mdata_(mesh,ic,jc,k);
            if (ic < width[BDRY2D_XHI])
                mesh->mdata_(mesh,ic+mesh->maxi, jc-mesh->maxj, k) = mesh->mdata_(mesh,ic,jc,k);
        }
    }
    
}
