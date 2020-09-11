/*************************************************
    quarterbc.cpp

    $Header: um2periodic.cpp,v 2.3 99/01/06 13:59:04 chopp Exp $

    $Log:       um2periodic.cpp,v $
Revision 2.3  99/01/06  13:59:04  13:59:04  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:05:49  10:05:49  chopp (David Chopp)
Initial revision

*************************************************/

#include "quarterbc.h"

namespace levelset {

    void UM2_QuarterBdry::Apply(const int k)
    {
        int i, j;
        
        for (i=-width[BDRY2D_XLO]; i<0; ++i)
            for (j=0; j<mesh->maxj; ++j) 
                mesh->mdata_(mesh,i,j,k) = mesh->mdata_(mesh,j,-i,k);

        for (i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i)
            for (j=0; j<mesh->maxj; ++j) 
                mesh->mdata_(mesh,i,j,k) = (i-mesh->maxi+2)*mesh->mdata_(mesh,mesh->maxi-1,j,k)
                    -(i-mesh->maxi+1)*mesh->mdata_(mesh,mesh->maxi-2,j,k);
   
        for (j=-width[BDRY2D_YLO]; j<0; ++j)
            for (i=0; i<mesh->maxi; ++i) 
                mesh->mdata_(mesh,i,j,k) = mesh->mdata_(mesh,-j,i,k);

        for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
            for (i=0; i<mesh->maxi; ++i) 
                mesh->mdata_(mesh,i,j,k) = (j-mesh->maxj+2)*mesh->mdata_(mesh,i,mesh->maxj-1,k)
                    -(j-mesh->maxj+1)*mesh->mdata_(mesh,i,mesh->maxj-2,k);

        for (i=-width[BDRY2D_XLO]; i<0; ++i) {
            for (j=-width[BDRY2D_YLO]; j<0; ++j)
                mesh->mdata_(mesh,i,j,k) = mesh->mdata_(mesh,-j,-i,k);
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
                mesh->mdata_(mesh,i,j,k) = mesh->mdata_(mesh,j,-i,k);
        }

        for (i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i) {
            for (j=-width[BDRY2D_YLO]; j<0; ++j)
                mesh->mdata_(mesh,i,j,k) = mesh->mdata_(mesh,-j,i,k);
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
                mesh->mdata_(mesh,i,j,k) = (j-mesh->maxj+1)/double(j-mesh->maxj+1-i)
                    *((i-mesh->maxi+2)*mesh->mdata_(mesh,mesh->maxi-1,j,k)
                      -(i-mesh->maxi+1)*mesh->mdata_(mesh,mesh->maxi-2,j,k))
                    -i/double(j-mesh->maxj+1-i)
                    *((j-mesh->maxj+2)*mesh->mdata_(mesh,i,mesh->maxj-1,k)
                      -(j-mesh->maxj+1)*mesh->mdata_(mesh,i,mesh->maxj-2,k));
        }
    }

    void UM2_QuarterBdry::Apply(const int ic, const int jc, const int k)
    {
        std::cerr << "Error: QuarterBdry individual apply not yet implemented\n";
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

}
