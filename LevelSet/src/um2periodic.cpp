/*************************************************
    um2periodic.cpp

    $Header: um2periodic.cpp,v 2.3 99/01/06 13:59:04 chopp Exp $

    $Log:       um2periodic.cpp,v $
Revision 2.3  99/01/06  13:59:04  13:59:04  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:05:49  10:05:49  chopp (David Chopp)
Initial revision

*************************************************/

#include "um2periodic.h"

namespace levelset {
        
    void UM2_PeriodicBdry::Apply(const int k)
    {
        int i, j;
   
        int xdiff = mesh->Index(mesh->maxi,0,0) - mesh->Index(0,0,0);
        int ydiff = mesh->Index(0,mesh->maxj,0) - mesh->Index(0,0,0);
   
        for (i=0; i<width[BDRY2D_XLO]; ++i)
            for (j=0; j<mesh->tmaxj; ++j)
                mesh->thedata[mesh->Index(i,j,k)]
                    = mesh->thedata[mesh->Index(i,j,k)+xdiff];

        for (i=0; i<width[BDRY2D_XHI]; ++i)
            for (j=0; j<mesh->tmaxj; ++j)
                mesh->thedata[mesh->Index(mesh->tmaxi-1-i,j,k)]
                    = mesh->thedata[mesh->Index(mesh->tmaxi-1-i,j,k)-xdiff];

        for (j=0; j<width[BDRY2D_YLO]; ++j)
            for (i=0; i<mesh->tmaxi; ++i)
                mesh->thedata[mesh->Index(i,j,k)]
                    = mesh->thedata[mesh->Index(i,j,k)+ydiff];

        for (j=0; j<width[BDRY2D_YHI]; ++j)
            for (i=0; i<mesh->tmaxi; ++i)
                mesh->thedata[mesh->Index(i,mesh->tmaxj-1-j,k)]
                    = mesh->thedata[mesh->Index(i,mesh->tmaxj-1-j,k)-ydiff];
    }

    void UM2_PeriodicBdry::Apply(const int ic, const int jc, const int k)
    {
        if (ic < width[BDRY2D_XHI]) {
            mesh->mdata_(mesh,ic+mesh->maxi, jc, k) = mesh->mdata_(mesh,ic,jc,k);
            if (jc < width[BDRY2D_YHI])
                mesh->mdata_(mesh,ic+mesh->maxi, jc+mesh->maxj, k) = mesh->mdata_(mesh,ic,jc,k);
        }
        if (jc < width[BDRY2D_YHI]) {
            mesh->mdata_(mesh,ic, jc+mesh->maxj, k) = mesh->mdata_(mesh,ic,jc,k);
            if (ic > mesh->maxi-1-width[BDRY2D_XLO])
                mesh->mdata_(mesh,ic-mesh->maxi, jc+mesh->maxj, k) = mesh->mdata_(mesh,ic,jc,k);
        }
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

    int UM2_PeriodicBdry::Nabor(const int i, const int j, const BCDir d, int& ni, int& nj)
    {
        int valid = 1;
        switch(d) {
        case Right: if (i < mesh->maxi) {ni=i+1; nj=j;} else {ni=0; nj=j;} break;
        case Left: if (i > 0) {ni=i-1; nj=j;} else {ni=mesh->maxi-1; nj=j;} break;
        case Up: if (j < mesh->maxj) {ni=i; nj=j+1;} else {ni=i; nj=0;} break;
        case Down: if (j > 0) {ni=i; nj=j-1;} else {ni=i; nj=mesh->maxj-1;}; break;
        }
        return valid;
    }

}
