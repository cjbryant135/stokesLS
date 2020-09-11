#include "um3boundary.h"

namespace levelset {
        
    void UM3_Boundary::SetIBoundary(const int l, const int n)
    {
        int i, j, k;

        for (i=-width[BDRY3D_XLO]; i<mesh->maxi+width[BDRY3D_XHI]; ++i) {
            for (j=-width[BDRY3D_YLO]; j<mesh->maxj+width[BDRY3D_YHI]; ++j) {
                for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                    mesh->midata_(mesh,i,j,k,l) = n;
                for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                    mesh->midata_(mesh,i,j,k,l) = n;
            }
            for (k=0; k<mesh->maxk; ++k) {
                for (j=-width[BDRY3D_YLO]; j<0; ++j)
                    mesh->midata_(mesh,i,j,k,l) = n;
                for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                    mesh->midata_(mesh,i,j,k,l) = n;
            }
        }
        for (i=-width[BDRY3D_XLO]; i<0; ++i)
            for (j=0; j<mesh->maxj; ++j)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->midata_(mesh,i,j,k,l) = n;
        for (i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i)
            for (j=0; j<mesh->maxj; ++j)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->midata_(mesh,i,j,k,l) = n;
    }
        
}
