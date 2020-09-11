/*************************************************
    um2periodic.cpp

    $Header: um2periodic.cpp,v 2.3 99/01/06 13:59:04 chopp Exp $

    $Log:       um2periodic.cpp,v $
Revision 2.3  99/01/06  13:59:04  13:59:04  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:05:49  10:05:49  chopp (David Chopp)
Initial revision

*************************************************/

#include "um2xperiodic.h"

namespace levelset {
        
    void UM2_XPeriodicBdry::Apply(const int k)
    {
        int i, j;
  		   
        for (j=-width[BDRY2D_YLO]; j<0; ++j) 
            for (i=0; i<mesh->maxi; ++i) 
                mesh->mdata_(mesh,i,j,k) = (1-j)*mesh->mdata_(mesh,i,0,k)
                    +j*mesh->mdata_(mesh,i,1,k);

        for (j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
            for (i=0; i<mesh->maxi; ++i) 
                mesh->mdata_(mesh,i,j,k) = (j-mesh->maxj+2)*mesh->mdata_(mesh,i,mesh->maxj-1,k)
                    -(j-mesh->maxj+1)*mesh->mdata_(mesh,i,mesh->maxj-2,k);

//OLD WAY
//        int xdiff = mesh->Index(mesh->maxi,0,0) - mesh->Index(0,0,0);
//   
//        for (i=0; i<width[BDRY2D_XLO]; ++i)
//            for (j=0; j<mesh->tmaxj; ++j)
//                mesh->thedata[mesh->Index(i,j,k)]
//                    = mesh->thedata[mesh->Index(i,j,k)+xdiff];
//
//        for (i=0; i<width[BDRY2D_XHI]; ++i)
//            for (j=0; j<mesh->tmaxj; ++j)
//                mesh->thedata[mesh->Index(mesh->tmaxi-1-i,j,k)]
//                    = mesh->thedata[mesh->Index(mesh->tmaxi-1-i,j,k)-xdiff];
			//NEW WAY (edited by Colton Bryant, Jun 23, 2020)
			for(i=-width[BDRY2D_XLO]; i<0; ++i)
				for(j=-width[BDRY2D_YLO];j<mesh->maxj+width[BDRY2D_YHI]; ++j)
					mesh->mdata_(mesh,i,j,k) = mesh->mdata_(mesh,mesh->maxi+i-1,j,k);
			
			for(i=mesh->maxi; i<mesh->maxi+width[BDRY2D_XHI]; ++i) 
				for(j=-width[BDRY2D_YLO];j<mesh->maxj+width[BDRY2D_YHI]; ++j)
					mesh->mdata_(mesh,i,j,k) = mesh->mdata_(mesh,i-mesh->maxi+1,j,k);
				

    }

    void UM2_XPeriodicBdry::Apply(const int ic, const int jc, const int k)
    {
        if (ic < width[BDRY2D_XHI]) {
            mesh->mdata_(mesh,ic+mesh->maxi, jc, k) = mesh->mdata_(mesh,ic,jc,k);
            if (jc < width[BDRY2D_YHI])
                mesh->mdata_(mesh,ic+mesh->maxi, jc+mesh->maxj, k) = mesh->mdata_(mesh,ic,jc,k);
        }
        if (ic > mesh->maxi-1-width[BDRY2D_XLO]) {
            mesh->mdata_(mesh,ic-mesh->maxi, jc, k) = mesh->mdata_(mesh,ic,jc,k);
            if (jc > mesh->maxj-1-width[BDRY2D_YLO])
                mesh->mdata_(mesh,ic-mesh->maxi, jc-mesh->maxj, k) = mesh->mdata_(mesh,ic,jc,k);
        }
        if (jc == 0 || jc == 1) 
            for (int j=-width[BDRY2D_YLO]; j < 0; ++j)
                mesh->mdata_(mesh,ic, j, k) = (1-j)*mesh->mdata_(mesh,ic,0,k)
                    +j*mesh->mdata_(mesh,ic,1,k);
        if (jc == mesh->maxj-1 || jc == mesh->maxj-2) 
            for (int j=mesh->maxj; j<mesh->maxj+width[BDRY2D_YHI]; ++j)
                mesh->mdata_(mesh,ic,j,k) = (j-mesh->maxj+2)*mesh->mdata_(mesh,ic,mesh->maxj-1,k)
                    -(j-mesh->maxj+1)*mesh->mdata_(mesh,ic,mesh->maxj-2,k);
    }

    int UM2_XPeriodicBdry::Nabor(const int i, const int j, const BCDir d, int& ni, int& nj)
    {
        int valid = 1;
        switch(d) {
        case Right: if (i < mesh->maxi) {ni=i+1; nj=j;} else {ni=0; nj=j;} break;
        case Left: if (i > 0) {ni=i-1; nj=j;} else {ni=mesh->maxi-1; nj=j;} break;
        case Up: if (j < mesh->maxj) {ni=i; nj=j+1;} else valid=0; break;
        case Down: if (j > 0) {ni=i; nj=j-1;} else valid=0; break;
        }
        return valid;
    }

}
