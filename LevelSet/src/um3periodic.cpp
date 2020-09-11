/*************************************************
    um3periodic.cpp

    $Header: um3periodic.cpp,v 2.3 99/01/06 13:59:09 chopp Exp $

    $Log:       um3periodic.cpp,v $
Revision 2.3  99/01/06  13:59:09  13:59:09  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:05:49  10:05:49  chopp (David Chopp)
Initial revision

*************************************************/

#include "um3periodic.h"

namespace levelset {
        
    void UM3_PeriodicBdry::Apply(const int l)
    {
        int i, j, k;
   
        int xdiff = mesh->Index(mesh->maxi,0,0,0) - mesh->Index(0,0,0,0);
        int ydiff = mesh->Index(0,mesh->maxj,0,0) - mesh->Index(0,0,0,0);
        int zdiff = mesh->Index(0,0,mesh->maxk,0) - mesh->Index(0,0,0,0);
   
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

        for (k=0; k<width[BDRY3D_ZLO]; ++k)
            for (i=0; i<mesh->tmaxi; ++i)
                for (j=0; j<mesh->tmaxj; ++j)
                    mesh->thedata[mesh->Index(i,j,k,l)]
                        = mesh->thedata[mesh->Index(i,j,k,l)+zdiff];

        for (k=0; k<width[BDRY3D_ZHI]; ++k)
            for (i=0; i<mesh->tmaxi; ++i)
                for (j=0; j<mesh->tmaxj; ++j)
                    mesh->thedata[mesh->Index(i,j,mesh->tmaxk-1-k,l)]
                        = mesh->thedata[mesh->Index(i,j,mesh->tmaxk-1-k,l)-zdiff];
    }

    void UM3_PeriodicBdry::Apply(const int ic, const int jc, const int kc,
                                 const int l)
    {
        if (ic < width[BDRY3D_XHI]) {
            mesh->data_(ic+mesh->maxi, jc, kc, l) = mesh->data_(ic,jc,kc,l);
            if (jc < width[BDRY3D_YHI]) {
                mesh->data_(ic+mesh->maxi, jc+mesh->maxj, kc, l)
                    = mesh->data_(ic, jc, kc, l);
                if (kc < width[BDRY3D_ZHI]) 
                    mesh->data_(ic+mesh->maxi, jc+mesh->maxj, kc+mesh->maxk, l)
                        = mesh->data_(ic, jc, kc, l);
                if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
                    mesh->data_(ic+mesh->maxi, jc+mesh->maxj, kc-mesh->maxk, l)
                        = mesh->data_(ic, jc, kc, l);
            }
            if (jc > mesh->maxj-1-width[BDRY3D_YLO]) {         
                mesh->data_(ic+mesh->maxi, jc-mesh->maxj, kc, l)
                    = mesh->data_(ic, jc, kc, l);
                if (kc < width[BDRY3D_ZHI]) 
                    mesh->data_(ic+mesh->maxi, jc-mesh->maxj, kc+mesh->maxk, l)
                        = mesh->data_(ic, jc, kc, l);
                if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
                    mesh->data_(ic+mesh->maxi, jc-mesh->maxj, kc-mesh->maxk, l)
                        = mesh->data_(ic, jc, kc, l);
            }
            if (kc < width[BDRY3D_ZHI]) 
                mesh->data_(ic+mesh->maxi, jc, kc+mesh->maxk, l)
                    = mesh->data_(ic, jc, kc, l);
            if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
                mesh->data_(ic+mesh->maxi, jc, kc-mesh->maxk, l)
                    = mesh->data_(ic, jc, kc, l);
        }
   
        if (ic > mesh->maxi-1-width[BDRY3D_XLO]) {
            mesh->data_(ic-mesh->maxi, jc, kc, l) = mesh->data_(ic,jc,kc,l);
            if (jc < width[BDRY3D_YHI]) {
                mesh->data_(ic-mesh->maxi, jc+mesh->maxj, kc, l)
                    = mesh->data_(ic, jc, kc, l);
                if (kc < width[BDRY3D_ZHI]) 
                    mesh->data_(ic-mesh->maxi, jc+mesh->maxj, kc+mesh->maxk, l)
                        = mesh->data_(ic, jc, kc, l);
                if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
                    mesh->data_(ic-mesh->maxi, jc+mesh->maxj, kc-mesh->maxk, l)
                        = mesh->data_(ic, jc, kc, l);
            }
            if (jc > mesh->maxj-1-width[BDRY3D_YLO]) {         
                mesh->data_(ic-mesh->maxi, jc-mesh->maxj, kc, l)
                    = mesh->data_(ic, jc, kc, l);
                if (kc < width[BDRY3D_ZHI]) 
                    mesh->data_(ic-mesh->maxi, jc-mesh->maxj, kc+mesh->maxk, l)
                        = mesh->data_(ic, jc, kc, l);
                if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
                    mesh->data_(ic-mesh->maxi, jc-mesh->maxj, kc-mesh->maxk, l)
                        = mesh->data_(ic, jc, kc, l);
            }
            if (kc < width[BDRY3D_ZHI]) 
                mesh->data_(ic-mesh->maxi, jc, kc+mesh->maxk, l)
                    = mesh->data_(ic, jc, kc, l);
            if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
                mesh->data_(ic-mesh->maxi, jc, kc-mesh->maxk, l)
                    = mesh->data_(ic, jc, kc, l);
        }

        if (jc < width[BDRY3D_YHI]) {
            mesh->data_(ic, jc+mesh->maxj, kc, l) = mesh->data_(ic, jc, kc, l);
            if (kc < width[BDRY3D_ZHI])
                mesh->data_(ic, jc+mesh->maxj, kc+mesh->maxk, l)
                    = mesh->data_(ic, jc, kc, l);
            if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
                mesh->data_(ic, jc+mesh->maxj, kc-mesh->maxk, l)
                    = mesh->data_(ic, jc, kc, l);
        }
   
        if (jc > mesh->maxj-1-width[BDRY3D_YLO]) {
            mesh->data_(ic, jc-mesh->maxj, kc, l) = mesh->data_(ic, jc, kc, l);
            if (kc < width[BDRY3D_ZHI])
                mesh->data_(ic, jc-mesh->maxj, kc+mesh->maxk, l)
                    = mesh->data_(ic, jc, kc, l);
            if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
                mesh->data_(ic, jc-mesh->maxj, kc-mesh->maxk, l)
                    = mesh->data_(ic, jc, kc, l);
        }

        if (kc < width[BDRY3D_ZHI]) 
            mesh->data_(ic, jc, kc+mesh->maxk, l) = mesh->data_(ic, jc, kc, l);
        if (kc > mesh->maxk-1-width[BDRY3D_ZLO])
            mesh->data_(ic, jc, kc-mesh->maxk, l) = mesh->data_(ic, jc, kc, l);
    }

}


