/*************************************************
    um3linear.cpp

    $Header: um3linear.cpp,v 2.3 99/01/06 13:59:08 chopp Exp $

    $Log:       um3linear.cpp,v $
Revision 2.3  99/01/06  13:59:08  13:59:08  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:06:16  10:06:16  chopp (David Chopp)
Initial revision

*************************************************/

#include "um3linear.h"

namespace levelset {
        
    void UM3_LinearBdry::Apply(const int l)
    {
        int i, j, k;

        for (i=-width[BDRY3D_XLO]; i<0; ++i)
            for (j=0; j<mesh->maxj; ++j)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->data_(i,j,k,l) = (1-i)*mesh->data_(0,j,k,l)
                        +i*mesh->data_(1,j,k,l);

        for (i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i)
            for (j=0; j<mesh->maxj; ++j)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->data_(i,j,k,l)
                        = (i-mesh->maxi+2)*mesh->data_(mesh->maxi-1,j,k,l)
                        -(i-mesh->maxi+1)*mesh->data_(mesh->maxi-2,j,k,l);
   
        for (j=-width[BDRY3D_YLO]; j<0; ++j)
            for (i=0; i<mesh->maxi; ++i)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->data_(i,j,k,l) = (1-j)*mesh->data_(i,0,k,l)
                        +j*mesh->data_(i,1,k,l);

        for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
            for (i=0; i<mesh->maxi; ++i)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->data_(i,j,k,l)
                        = (j-mesh->maxj+2)*mesh->data_(i,mesh->maxj-1,k,l)
                        -(j-mesh->maxj+1)*mesh->data_(i,mesh->maxj-2,k,l);

        for (k=-width[BDRY3D_ZLO]; k<0; ++k)
            for (i=0; i<mesh->maxi; ++i)
                for (j=0; j<mesh->maxj; ++j)
                    mesh->data_(i,j,k,l) = (1-k)*mesh->data_(i,j,0,l)
                        +k*mesh->data_(i,j,1,l);

        for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
            for (i=0; i<mesh->maxi; ++i)
                for (j=0; j<mesh->maxj; ++j)
                    mesh->data_(i,j,k,l)
                        = (k-mesh->maxk+2)*mesh->data_(i,j,mesh->maxk-1,l)
                        -(k-mesh->maxk+1)*mesh->data_(i,j,mesh->maxk-2,l);

        for (i=-width[BDRY3D_XLO]; i<0; ++i) {
            for (j=-width[BDRY3D_YLO]; j<0; ++j)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->data_(i,j,k,l) = i*mesh->data_(1,0,k,l)
                        +j*mesh->data_(0,1,k,l)+(1-i-j)*mesh->data_(0,0,k,l);
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->data_(i,j,k,l)
                        = i*mesh->data_(1,mesh->maxj-1,k,l)
                        +(mesh->maxj-1-j)*mesh->data_(0,mesh->maxj-2,k,l)
                        +(2-mesh->maxj-i+j)*mesh->data_(0,mesh->maxj-1,k,l);
            for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                for (j=0; j<mesh->maxj; ++j)
                    mesh->data_(i,j,k,l) = i*mesh->data_(1,j,0,l)
                        +k*mesh->data_(0,j,1,l)+(1-i-k)*mesh->data_(0,j,0,l);
            for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                for (j=0; j<mesh->maxj; ++j)
                    mesh->data_(i,j,k,l) 
                        = i*mesh->data_(1,j,mesh->maxk-1,l)
                        +(mesh->maxk-1-k)*mesh->data_(0,j,mesh->maxk-2,l)
                        +(2-mesh->maxk-i+k)*mesh->data_(0,j,mesh->maxk-1,l);
      
            for (j=-width[BDRY3D_YLO]; j<0; ++j)
                for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                    mesh->data_(i,j,k,l) = i*mesh->data_(1,0,0,l)
                        +j*mesh->data_(0,1,0,l)+k*mesh->data_(0,0,1,l)
                        +(1-i-j-k)*mesh->data_(0,0,0,l);
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                    mesh->data_(i,j,k,l) = i*mesh->data_(1,mesh->maxj-1,0,l)
                        +(mesh->maxj-1-j)*mesh->data_(0,mesh->maxj-2,0,l)
                        +k*mesh->data_(0,mesh->maxj-1,1,l)
                        +(2-mesh->maxj-i+j-k)*mesh->data_(0,mesh->maxj-1,0,l);
            for (j=-width[BDRY3D_YLO]; j<0; ++j)
                for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                    mesh->data_(i,j,k,l) = i*mesh->data_(1,0,mesh->maxk-1,l)
                        +j*mesh->data_(0,1,mesh->maxk-1,l)
                        +(mesh->maxk-1-k)*mesh->data_(0,0,mesh->maxk-2,l)
                        +(2-mesh->maxk-i-j+k)*mesh->data_(0,0,mesh->maxk-1,l);
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                    mesh->data_(i,j,k,l) = i*mesh->data_(1,mesh->maxj-1,mesh->maxk-1,l)
                        +(mesh->maxj-1-j)*mesh->data_(0,mesh->maxj-2,mesh->maxk-1,l)
                        +(mesh->maxk-1-k)*mesh->data_(0,mesh->maxj-1,mesh->maxk-2,l)
                        +(3-mesh->maxj-mesh->maxk-i+j+k)
                        *mesh->data_(0,mesh->maxj-1,mesh->maxk-1,l);
        }
        for (i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) {
            for (j=-width[BDRY3D_YLO]; j<0; ++j)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->data_(i,j,k,l)
                        = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,0,k,l)
                        +j*mesh->data_(mesh->maxi-1,1,k,l)
                        +(2-mesh->maxi+i-j)*mesh->data_(mesh->maxi-1,0,k,l);
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                for (k=0; k<mesh->maxk; ++k)
                    mesh->data_(i,j,k,l)
                        = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,mesh->maxj-1,k,l)
                        +(mesh->maxj-1-j)*mesh->data_(mesh->maxi-1,mesh->maxj-2,k,l)
                        +(3-mesh->maxj-mesh->maxi+i+j)
                        *mesh->data_(mesh->maxi-1,mesh->maxj-1,k,l);
            for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                for (j=0; j<mesh->maxj; ++j)
                    mesh->data_(i,j,k,l)
                        = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,j,0,l)
                        +k*mesh->data_(mesh->maxi-1,j,1,l)
                        +(2-mesh->maxi+i-k)*mesh->data_(mesh->maxi-1,j,0,l);
            for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                for (j=0; j<mesh->maxj; ++j)
                    mesh->data_(i,j,k,l) 
                        = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,j,mesh->maxk-1,l)
                        +(mesh->maxk-1-k)*mesh->data_(mesh->maxi-1,j,mesh->maxk-2,l)
                        +(3-mesh->maxk-mesh->maxi+i+k)
                        *mesh->data_(mesh->maxi-1,j,mesh->maxk-1,l);
      
            for (j=-width[BDRY3D_YLO]; j<0; ++j)
                for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                    mesh->data_(i,j,k,l)
                        = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,0,0,l)
                        +j*mesh->data_(mesh->maxi-1,1,0,l)
                        +k*mesh->data_(mesh->maxi-1,0,1,l)
                        +(2-mesh->maxi+i-j-k)*mesh->data_(mesh->maxi-1,0,0,l);
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                    mesh->data_(i,j,k,l)
                        = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,mesh->maxj-1,0,l)
                        +(mesh->maxj-1-j)*mesh->data_(mesh->maxi-1,mesh->maxj-2,0,l)
                        +k*mesh->data_(mesh->maxi-1,mesh->maxj-1,1,l)
                        +(3-mesh->maxj-mesh->maxi+i+j-k)
                        *mesh->data_(mesh->maxi-1,mesh->maxj-1,0,l);
            for (j=-width[BDRY3D_YLO]; j<0; ++j)
                for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                    mesh->data_(i,j,k,l)
                        = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,0,mesh->maxk-1,l)
                        +j*mesh->data_(mesh->maxi-1,1,mesh->maxk-1,l)
                        +(mesh->maxk-1-k)*mesh->data_(mesh->maxi-1,0,mesh->maxk-2,l)
                        +(3-mesh->maxk-mesh->maxi+i-j+k)
                        *mesh->data_(mesh->maxi-1,0,mesh->maxk-1,l);
            for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                    mesh->data_(i,j,k,l)
                        = (mesh->maxi-1-i)
                        *mesh->data_(mesh->maxi-2,mesh->maxj-1,mesh->maxk-1,l)
                        +(mesh->maxj-1-j)
                        *mesh->data_(mesh->maxi-1,mesh->maxj-2,mesh->maxk-1,l)
                        +(mesh->maxk-1-k)
                        *mesh->data_(mesh->maxi-1,mesh->maxj-1,mesh->maxk-2,l)
                        +(4-mesh->maxj-mesh->maxk-mesh->maxi+i+j+k)
                        *mesh->data_(mesh->maxi-1,mesh->maxj-1,mesh->maxk-1,l);
        }
      
        for (j=-width[BDRY3D_YLO]; j<0; ++j) {
            for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                for (i=0; i<mesh->maxi; ++i)
                    mesh->data_(i,j,k,l) = j*mesh->data_(i,1,0,l)
                        +k*mesh->data_(i,0,1,l)+(1-j-k)*mesh->data_(i,0,0,l);
            for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                for (i=0; i<mesh->maxi; ++i)
                    mesh->data_(i,j,k,l)
                        = j*mesh->data_(i,1,mesh->maxk-1,l)
                        +(mesh->maxk-1-k)*mesh->data_(i,0,mesh->maxk-2,l)
                        +(2-mesh->maxk-j+k)*mesh->data_(i,0,mesh->maxk-1,l);
        }
        for (j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j) {
            for (k=-width[BDRY3D_ZLO]; k<0; ++k)
                for (i=0; i<mesh->maxi; ++i)
                    mesh->data_(i,j,k,l)
                        = (mesh->maxj-1-j)*mesh->data_(i,mesh->maxi-2,0,l)
                        +k*mesh->data_(i,mesh->maxj-1,1,l)
                        +(2-mesh->maxj+j-k)*mesh->data_(i,mesh->maxj-1,0,l);
            for (k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                for (i=0; i<mesh->maxi; ++i)
                    mesh->data_(i,j,k,l) 
                        = (mesh->maxj-1-j)*mesh->data_(i,mesh->maxj-2,mesh->maxk-1,l)
                        +(mesh->maxk-1-k)*mesh->data_(i,mesh->maxj-1,mesh->maxk-2,l)
                        +(3-mesh->maxk-mesh->maxj+j+k)
                        *mesh->data_(i,mesh->maxj-1,mesh->maxk-1,l);
        }
      
    }

    void UM3_LinearBdry::Apply(const int ic, const int jc, const int kc,
                               const int l)
    {
        if (ic == 0 || ic == 1) {
            for (int i=-width[BDRY3D_XLO]; i<0; ++i)
                mesh->data_(i,jc,kc,l) = (1-i)*mesh->data_(0,jc,kc,l)
                    +i*mesh->data_(1,jc,kc,l);
            if (jc == 0 || jc == 1) {
                for (int i=-width[BDRY3D_XLO]; i<0; ++i) 
                    for (int j=-width[BDRY3D_YLO]; j<0; ++j)
                        mesh->data_(i,j,kc,l) = i*mesh->data_(1,0,kc,l)
                            +j*mesh->data_(0,1,kc,l)+(1-i-j)*mesh->data_(0,0,kc,l);
                if (kc == 0 || kc == 1) 
                    for (int i=-width[BDRY3D_XLO]; i<0; ++i) 
                        for (int j=-width[BDRY3D_YLO]; j<0; ++j)
                            for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                                mesh->data_(i,j,k,l) = i*mesh->data_(1,0,0,l)
                                    +j*mesh->data_(0,1,0,l)+k*mesh->data_(0,0,1,l)
                                    +(1-i-j-k)*mesh->data_(0,0,0,l);
                if (kc == mesh->maxk-1 || kc == mesh->maxk-2)
                    for (int i=-width[BDRY3D_XLO]; i<0; ++i) 
                        for (int j=-width[BDRY3D_YLO]; j<0; ++j)
                            for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                                mesh->data_(i,j,k,l) = i*mesh->data_(1,0,mesh->maxk-1,l)
                                    +j*mesh->data_(0,1,mesh->maxk-1,l)
                                    +(mesh->maxk-1-k)*mesh->data_(0,0,mesh->maxk-2,l)
                                    +(2-mesh->maxk-i-j+k)
                                    *mesh->data_(0,0,mesh->maxk-1,l);
            }
            if (jc == mesh->maxj-1 || jc == mesh->maxj-2) {
                for (int i=-width[BDRY3D_XLO]; i<0; ++i) 
                    for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                        mesh->data_(i,j,kc,l)
                            = i*mesh->data_(1,mesh->maxj-1,kc,l)
                            +(mesh->maxj-1-j)*mesh->data_(0,mesh->maxj-2,kc,l)
                            +(2-mesh->maxj-i+j)*mesh->data_(0,mesh->maxj-1,kc,l);
                if (kc == 0 || kc == 1) 
                    for (int i=-width[BDRY3D_XLO]; i<0; ++i) 
                        for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                            for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                                mesh->data_(i,j,k,l) = i*mesh->data_(1,mesh->maxj-1,0,l)
                                    +(mesh->maxj-1-j)*mesh->data_(0,mesh->maxj-2,0,l)
                                    +k*mesh->data_(0,mesh->maxj-1,1,l)
                                    +(2-mesh->maxj-i+j-k)
                                    *mesh->data_(0,mesh->maxj-1,0,l);
                if (kc == mesh->maxk-1 || kc == mesh->maxk-2)
                    for (int i=-width[BDRY3D_XLO]; i<0; ++i) 
                        for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                            for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                                mesh->data_(i,j,k,l)
                                    = i*mesh->data_(1,mesh->maxj-1,mesh->maxk-1,l)
                                    +(mesh->maxj-1-j)
                                    *mesh->data_(0,mesh->maxj-2,mesh->maxk-1,l)
                                    +(mesh->maxk-1-k)
                                    *mesh->data_(0,mesh->maxj-1,mesh->maxk-2,l)
                                    +(3-mesh->maxj-mesh->maxk-i+j+k)
                                    *mesh->data_(0,mesh->maxj-1,mesh->maxk-1,l);
            }
            if (kc == 0 || kc == 1) 
                for (int i=-width[BDRY3D_XLO]; i<0; ++i) 
                    for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                        mesh->data_(i,jc,k,l) = i*mesh->data_(1,jc,0,l)
                            +k*mesh->data_(0,jc,1,l)+(1-i-k)*mesh->data_(0,jc,0,l);
      
            if (kc == mesh->maxk-1 || kc == mesh->maxk-2) 
                for (int i=-width[BDRY3D_XLO]; i<0; ++i) 
                    for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                        mesh->data_(i,jc,k,l) 
                            = i*mesh->data_(1,jc,mesh->maxk-1,l)
                            +(mesh->maxk-1-k)*mesh->data_(0,jc,mesh->maxk-2,l)
                            +(2-mesh->maxk-i+k)*mesh->data_(0,jc,mesh->maxk-1,l);
        }
   
        if (ic == mesh->maxi-1 || ic == mesh->maxi-2) {
            for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                mesh->data_(i,jc,kc,l)
                    = (i-mesh->maxi+2)*mesh->data_(mesh->maxi-1,jc,kc,l)
                    -(i-mesh->maxi+1)*mesh->data_(mesh->maxi-2,jc,kc,l);
            if (jc == 0 || jc == 1) {
                for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                    for (int j=-width[BDRY3D_YLO]; j<0; ++j)
                        mesh->data_(i,j,kc,l)
                            = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,0,kc,l)
                            +j*mesh->data_(mesh->maxi-1,1,kc,l)
                            +(2-mesh->maxi+i-j)*mesh->data_(mesh->maxi-1,0,kc,l);
                if (kc == 0 || kc == 1) 
                    for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                        for (int j=-width[BDRY3D_YLO]; j<0; ++j)
                            for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                                mesh->data_(i,j,k,l)
                                    = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,0,0,l)
                                    +j*mesh->data_(mesh->maxi-1,1,0,l)
                                    +k*mesh->data_(mesh->maxi-1,0,1,l)
                                    +(2-mesh->maxi+i-j-k)
                                    *mesh->data_(mesh->maxi-1,0,0,l);
                if (kc == mesh->maxk-1 || kc == mesh->maxk-2)
                    for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                        for (int j=-width[BDRY3D_YLO]; j<0; ++j)
                            for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                                mesh->data_(i,j,k,l)
                                    = (mesh->maxi-1-i)
                                    *mesh->data_(mesh->maxi-2,0,mesh->maxk-1,l)
                                    +j*mesh->data_(mesh->maxi-1,1,mesh->maxk-1,l)
                                    +(mesh->maxk-1-k)
                                    *mesh->data_(mesh->maxi-1,0,mesh->maxk-2,l)
                                    +(3-mesh->maxk-mesh->maxi+i-j+k)
                                    *mesh->data_(mesh->maxi-1,0,mesh->maxk-1,l);
            }
            if (jc == mesh->maxj-1 || jc == mesh->maxj-2) {
                for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                    for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                        mesh->data_(i,j,kc,l)
                            = (mesh->maxi-1-i)
                            *mesh->data_(mesh->maxi-2,mesh->maxj-1,kc,l)
                            +(mesh->maxj-1-j)
                            *mesh->data_(mesh->maxi-1,mesh->maxj-2,kc,l)
                            +(3-mesh->maxj-mesh->maxi+i+j)
                            *mesh->data_(mesh->maxi-1,mesh->maxj-1,kc,l);
                if (kc == 0 || kc == 1) 
                    for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                        for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                            for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                                mesh->data_(i,j,k,l)
                                    = (mesh->maxi-1-i)
                                    *mesh->data_(mesh->maxi-2,mesh->maxj-1,0,l)
                                    +(mesh->maxj-1-j)
                                    *mesh->data_(mesh->maxi-1,mesh->maxj-2,0,l)
                                    +k*mesh->data_(mesh->maxi-1,mesh->maxj-1,1,l)
                                    +(3-mesh->maxj-mesh->maxi+i+j-k)
                                    *mesh->data_(mesh->maxi-1,mesh->maxj-1,0,l);
                if (kc == mesh->maxk-1 || kc == mesh->maxk-2)
                    for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                        for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j)
                            for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                                mesh->data_(i,j,k,l)
                                    = (mesh->maxi-1-i)
                                    *mesh->data_(mesh->maxi-2,mesh->maxj-1,
                                                 mesh->maxk-1,l)
                                    +(mesh->maxj-1-j)
                                    *mesh->data_(mesh->maxi-1,mesh->maxj-2,
                                                 mesh->maxk-1,l)
                                    +(mesh->maxk-1-k)
                                    *mesh->data_(mesh->maxi-1,mesh->maxj-1,
                                                 mesh->maxk-2,l)
                                    +(4-mesh->maxj-mesh->maxk-mesh->maxi+i+j+k)
                                    *mesh->data_(mesh->maxi-1,mesh->maxj-1,
                                                 mesh->maxk-1,l);
            }
            if (kc == 0 || kc == 1) 
                for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                    for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                        mesh->data_(i,jc,k,l)
                            = (mesh->maxi-1-i)*mesh->data_(mesh->maxi-2,jc,0,l)
                            +k*mesh->data_(mesh->maxi-1,jc,1,l)
                            +(2-mesh->maxi+i-k)*mesh->data_(mesh->maxi-1,jc,0,l);
            if (kc == mesh->maxk-1 || kc == mesh->maxk-2) 
                for (int i=mesh->maxi; i<mesh->maxi+width[BDRY3D_XHI]; ++i) 
                    for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                        mesh->data_(i,jc,k,l) 
                            = (mesh->maxi-1-i)
                            *mesh->data_(mesh->maxi-2,jc,mesh->maxk-1,l)
                            +(mesh->maxk-1-k)
                            *mesh->data_(mesh->maxi-1,jc,mesh->maxk-2,l)
                            +(3-mesh->maxk-mesh->maxi+i+k)
                            *mesh->data_(mesh->maxi-1,jc,mesh->maxk-1,l);
        }

        if (jc == 0 || jc == 1) {
            for (int j=-width[BDRY3D_YLO]; j<0; ++j)
                mesh->data_(ic,j,kc,l) = (1-j)*mesh->data_(ic,0,kc,l)
                    +j*mesh->data_(ic,1,kc,l);

            if (kc == 0 || kc == 1) 
                for (int j=-width[BDRY3D_YLO]; j<0; ++j) 
                    for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                        mesh->data_(ic,j,k,l) = j*mesh->data_(ic,1,0,l)
                            +k*mesh->data_(ic,0,1,l)+(1-j-k)*mesh->data_(ic,0,0,l);

            if (kc == mesh->maxk-1 || kc == mesh->maxk-2)
                for (int j=-width[BDRY3D_YLO]; j<0; ++j) 
                    for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                        mesh->data_(ic,j,k,l)
                            = j*mesh->data_(ic,1,mesh->maxk-1,l)
                            +(mesh->maxk-1-k)*mesh->data_(ic,0,mesh->maxk-2,l)
                            +(2-mesh->maxk-j+k)*mesh->data_(ic,0,mesh->maxk-1,l);
        }
   
        if (jc == mesh->maxj-1 || jc == mesh->maxj-2) {
            for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j) 
                mesh->data_(ic,j,kc,l)
                    = (j-mesh->maxj+2)*mesh->data_(ic,mesh->maxj-1,kc,l)
                    -(j-mesh->maxj+1)*mesh->data_(ic,mesh->maxj-2,kc,l);

            if (kc == 0 || kc == 1) 
                for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j) 
                    for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                        mesh->data_(ic,j,k,l)
                            = (mesh->maxj-1-j)*mesh->data_(ic,mesh->maxj-2,0,l)
                            +k*mesh->data_(ic,mesh->maxj-1,1,l)
                            +(2-mesh->maxj+j-k)*mesh->data_(ic,mesh->maxj-1,0,l);
            if (kc == mesh->maxk-1 || kc == mesh->maxk-2)
                for (int j=mesh->maxj; j<mesh->maxj+width[BDRY3D_YHI]; ++j) 
                    for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                        mesh->data_(ic,j,k,l) 
                            = (mesh->maxj-1-j)
                            *mesh->data_(ic,mesh->maxj-2,mesh->maxk-1,l)
                            +(mesh->maxk-1-k)
                            *mesh->data_(ic,mesh->maxj-1,mesh->maxk-2,l)
                            +(3-mesh->maxk-mesh->maxj+j+k)
                            *mesh->data_(ic,mesh->maxj-1,mesh->maxk-1,l);
        }

        if (kc == 0 || kc == 1) 
            for (int k=-width[BDRY3D_ZLO]; k<0; ++k)
                mesh->data_(ic,jc,k,l) = (1-k)*mesh->data_(ic,jc,0,l)
                    +k*mesh->data_(ic,jc,1,l);
      
        if (kc == mesh->maxk-1 || kc == mesh->maxk-2) 
            for (int k=mesh->maxk; k<mesh->maxk+width[BDRY3D_ZHI]; ++k)
                mesh->data_(ic,jc,k,l)
                    = (k-mesh->maxk+2)*mesh->data_(ic,jc,mesh->maxk-1,l)
                    -(k-mesh->maxk+1)*mesh->data_(ic,jc,mesh->maxk-2,l);
    }

}
