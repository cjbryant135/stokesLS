/*************************************************
    um2linear.cpp

    $Header: um2linear.cpp,v 2.3 99/01/06 13:59:04 chopp Exp $

    $Log:       um2linear.cpp,v $
Revision 2.3  99/01/06  13:59:04  13:59:04  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:06:16  10:06:16  chopp (David Chopp)
Initial revision

*************************************************/

#include "um2linear.h"

namespace levelset {
        
    void UM2_LinearBdry::Apply(const int k)
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
   
        for (j=-width[BDRY2D_YLO]; j<0; ++j)
            for (i=0; i<mesh->maxi; ++i) 
                mesh->mdata_(mesh,i,j,k) = (1-j)*mesh->mdata_(mesh,i,0,k)
                    +j*mesh->mdata_(mesh,i,1,k);

        for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
            for (i=0; i<mesh->maxi; ++i) 
                mesh->mdata_(mesh,i,j,k) = (j-mesh->maxj+2)*mesh->mdata_(mesh,i,mesh->maxj-1,k)
                    -(j-mesh->maxj+1)*mesh->mdata_(mesh,i,mesh->maxj-2,k);

        for (i=-width[BDRY2D_XLO]; i<0; ++i) {
            for (j=-width[BDRY2D_YLO]; j<0; ++j)
                mesh->mdata_(mesh,i,j,k) = j/double(i+j)*((1-i)*mesh->mdata_(mesh,0,j,k)
                                                          +i*mesh->mdata_(mesh,1,j,k))
                    +i/double(i+j)*((1-j)*mesh->mdata_(mesh,i,0,k)
                                    +j*mesh->mdata_(mesh,i,1,k));
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
                mesh->mdata_(mesh,i,j,k) = (j-mesh->maxj+1)/double(j-mesh->maxj+1-i)
                    *((1-i)*mesh->mdata_(mesh,0,j,k)+i*mesh->mdata_(mesh,1,j,k))
                    -i/double(j-mesh->maxj+1-i)
                    *((j-mesh->maxj+2)*mesh->mdata_(mesh,i,mesh->maxj-1,k)
                      -(j-mesh->maxj+1)*mesh->mdata_(mesh,i,mesh->maxj-2,k));
        }

        for (i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i) {
            for (j=-width[BDRY2D_YLO]; j<0; ++j)
                mesh->mdata_(mesh,i,j,k) = j/double(i+j)
                    *((i-mesh->maxi+2)*mesh->mdata_(mesh,mesh->maxi-1,j,k)
                      -(i-mesh->maxi+1)*mesh->mdata_(mesh,mesh->maxi-2,j,k))
                    +i/double(i+j)*((1-j)*mesh->mdata_(mesh,i,0,k)
                                    +j*mesh->mdata_(mesh,i,1,k));
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
                mesh->mdata_(mesh,i,j,k) = (j-mesh->maxj+1)/double(j-mesh->maxj+1-i)
                    *((i-mesh->maxi+2)*mesh->mdata_(mesh,mesh->maxi-1,j,k)
                      -(i-mesh->maxi+1)*mesh->mdata_(mesh,mesh->maxi-2,j,k))
                    -i/double(j-mesh->maxj+1-i)
                    *((j-mesh->maxj+2)*mesh->mdata_(mesh,i,mesh->maxj-1,k)
                      -(j-mesh->maxj+1)*mesh->mdata_(mesh,i,mesh->maxj-2,k));
        }
    }

    void UM2_LinearBdry::Apply(const int ic, const int jc, const int k)
    {
        int i, j;
   
        if (ic == 0 || ic == 1) 
            for (i=-width[BDRY2D_XLO]; i < 0; ++i)
                mesh->mdata_(mesh,i, jc, k) = (1-i)*mesh->mdata_(mesh,0,jc,k)
                    +i*mesh->mdata_(mesh,1,jc,k);
        if (ic == mesh->maxi-1 || ic == mesh->maxi-2) 
            for (i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i)
                mesh->mdata_(mesh,i,jc,k) = (i-mesh->maxi+2)*mesh->mdata_(mesh,mesh->maxi-1,jc,k)
                    -(i-mesh->maxi+1)*mesh->mdata_(mesh,mesh->maxi-2,jc,k);
        if (jc == 0 || jc == 1) 
            for (j=-width[BDRY2D_YLO]; j < 0; ++j)
                mesh->mdata_(mesh,ic, j, k) = (1-j)*mesh->mdata_(mesh,ic,0,k)
                    +j*mesh->mdata_(mesh,ic,1,k);
        if (jc == mesh->maxj-1 || jc == mesh->maxj-2) 
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
                mesh->mdata_(mesh,ic,j,k) = (j-mesh->maxj+2)*mesh->mdata_(mesh,ic,mesh->maxj-1,k)
                    -(j-mesh->maxj+1)*mesh->mdata_(mesh,ic,mesh->maxj-2,k);
    }

}
