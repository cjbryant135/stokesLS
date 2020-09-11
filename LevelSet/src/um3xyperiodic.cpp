/*************************************************
    um3xyperiodic.cpp

    $Header: um3xyperiodic.cpp,v 2.3 99/01/06 13:59:09 chopp Exp $

    $Log:       um3xyperiodic.cpp,v $
Revision 2.3  99/01/06  13:59:09  13:59:09  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:05:49  10:05:49  chopp (David Chopp)
Initial revision

*************************************************/

#include "um3xyperiodic.h"

namespace levelset {
        
    void UM3_XYPeriodicBdry::Apply(const int l)
    {
        int i, j, k;
   
        int xdiff = mesh->Index(mesh->maxi,0,0,0) - mesh->Index(0,0,0,0);
        int ydiff = mesh->Index(0,mesh->maxj,0,0) - mesh->Index(0,0,0,0);
   
        for (i=0; i<width[BDRY3D_XLO]; ++i)
            for (j=0; j<mesh->tmaxj; ++j)
                for (k=0; k<mesh->tmaxk; ++k)
                    mesh->thedata[mesh->Index(i,j,k,l)]
                        = mesh->thedata[mesh->Index(i,j,k,l)+xdiff];

        for (i=0; i<width[BDRY3D_XHI]; ++i)
            for (j=0; j<mesh->tmaxj; ++j)
                for (k=0; k<mesh->tmaxk; ++k)
                    mesh->thedata[mesh->Index(mesh->tmaxi-1-i,j,k,l)]
                        = mesh->thedata[mesh->Index(mesh->tmaxi-1-i,j,k,l)-xdiff];

        for (j=0; j<width[BDRY3D_YLO]; ++j)
            for (i=0; i<mesh->tmaxi; ++i)
                for (k=0; k<mesh->tmaxk; ++k)
                    mesh->thedata[mesh->Index(i,j,k,l)]
                        = mesh->thedata[mesh->Index(i,j,k,l)+ydiff];

        for (j=0; j<width[BDRY3D_YHI]; ++j)
            for (i=0; i<mesh->tmaxi; ++i)
                for (k=0; k<mesh->tmaxk; ++k)
                    mesh->thedata[mesh->Index(i,mesh->tmaxj-1-j,k,l)]
                        = mesh->thedata[mesh->Index(i,mesh->tmaxj-1-j,k,l)-ydiff];

        for (k=-width[BDRY3D_ZLO]; k<0; ++k)
            for (i=-width[BDRY3D_XLO]; i<mesh->maxi+width[BDRY3D_XHI]; ++i)
                for (j=-width[BDRY3D_YLO]; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                    mesh->data_(i,j,k,l) = (1-k)*mesh->data_(i,j,0,l)
                        +k*mesh->data_(i,j,1,l);

        for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
            for (i=-width[BDRY3D_XLO]; i<mesh->maxi+width[BDRY3D_XHI]; ++i)
                for (j=-width[BDRY3D_YLO]; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                    mesh->data_(i,j,k,l)
                        = (k-mesh->maxk+2)*mesh->data_(i,j,mesh->maxk-1,l)
                        -(k-mesh->maxk+1)*mesh->data_(i,j,mesh->maxk-2,l);

    }

    void UM3_XYPeriodicBdry::Apply(const int ic, const int jc, const int kc,
                                   const int l)
    {
        if (kc < 2) {
            for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                mesh->data_(ic,jc,k,l) = (1-k)*mesh->data_(ic,jc,0,l)
                    +k*mesh->data_(ic,jc,1,l);
        }
        if (kc >= mesh->maxi-2) {
            for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                mesh->data_(ic,jc,k,l) = (k-mesh->maxk+2)*mesh->data_(ic,jc,mesh->maxk-1,l)
                    -(k-mesh->maxk+1)*mesh->data_(ic,jc,mesh->maxk-2,l);
        }
                
        if (ic < width[BDRY3D_XHI]) {
            mesh->data_(ic+mesh->maxi, jc, kc, l) = mesh->data_(ic,jc,kc,l);
            if (jc < width[BDRY3D_YHI]) {
                mesh->data_(ic+mesh->maxi, jc+mesh->maxj, kc, l)
                    = mesh->data_(ic, jc, kc, l);
            }
            if (jc > mesh->maxj-1-width[BDRY3D_YLO]) {         
                mesh->data_(ic+mesh->maxi, jc-mesh->maxj, kc, l)
                    = mesh->data_(ic, jc, kc, l);
            }
        }
   
        if (ic > mesh->maxi-1-width[BDRY3D_XLO]) {
            mesh->data_(ic-mesh->maxi, jc, kc, l) = mesh->data_(ic,jc,kc,l);
            if (jc < width[BDRY3D_YHI]) {
                mesh->data_(ic-mesh->maxi, jc+mesh->maxj, kc, l)
                    = mesh->data_(ic, jc, kc, l);
            }
            if (jc > mesh->maxj-1-width[BDRY3D_YLO]) {         
                mesh->data_(ic-mesh->maxi, jc-mesh->maxj, kc, l)
                    = mesh->data_(ic, jc, kc, l);
            }
        }

        if (jc < width[BDRY3D_YHI]) {
            mesh->data_(ic, jc+mesh->maxj, kc, l) = mesh->data_(ic, jc, kc, l);
        }
   
        if (jc > mesh->maxj-1-width[BDRY3D_YLO]) {
            mesh->data_(ic, jc-mesh->maxj, kc, l) = mesh->data_(ic, jc, kc, l);
        }

    }

}


