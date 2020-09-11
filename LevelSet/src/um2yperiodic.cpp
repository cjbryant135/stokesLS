/*************************************************
    um2yperiodic.cpp

    $Header: um2periodic.cpp,v 2.3 99/01/06 13:59:04 chopp Exp $

    $Log:       um2periodic.cpp,v $
Revision 2.3  99/01/06  13:59:04  13:59:04  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:05:49  10:05:49  chopp (David Chopp)
Initial revision

*************************************************/

#include "um2yperiodic.h"

namespace levelset {
        
    void UM2_YPeriodicBdry::Apply(const int k)
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
                mesh->thedata[mesh->Index(i,j,k)]
                    = mesh->thedata[mesh->Index(i,j,k)+ydiff];

        for (j=0; j<width[BDRY2D_YHI]; ++j)
            for (i=0; i<mesh->tmaxi; ++i)
                mesh->thedata[mesh->Index(i,mesh->tmaxj-1-j,k)]
                    = mesh->thedata[mesh->Index(i,mesh->tmaxj-1-j,k)-ydiff];
    }

    void UM2_YPeriodicBdry::Apply(const int ic, const int jc, const int k)
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

    int UM2_YPeriodicBdry::Nabor(const int i, const int j, const BCDir d, int& ni, int& nj)
    {
        int valid = 1;
        switch(d) {
        case Right: if (i < mesh->maxi) {ni=i+1; nj=j;} else valid=0; break;
        case Left: if (i > 0) {ni=i-1; nj=j;} else valid=0; break;
        case Up: if (j < mesh->maxj) {ni=i; nj=j+1;} else {ni=i; nj=0;} break;
        case Down: if (j > 0) {ni=i; nj=j-1;} else {ni=i; nj=mesh->maxj-1;}; break;
        }
        return valid;
    }

}
