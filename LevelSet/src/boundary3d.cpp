/*************************************************
    boundary3d.cpp

    $Header: boundary3d.cpp,v 2.3 99/01/06 13:59:05 chopp Exp $

    $Log:       boundary3d.cpp,v $
Revision 2.3  99/01/06  13:59:05  13:59:05  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:05:20  10:05:20  chopp (David Chopp)
Initial revision

*************************************************/

#include "boundary3d.h"

namespace levelset {
        
    Boundary3D::Boundary3D(const int wxlo, const int wxhi,
                           const int wylo, const int wyhi,
                           const int wzlo, const int wzhi)
    {
        width[BDRY3D_XLO] = wxlo;
        width[BDRY3D_XHI] = wxhi;
        width[BDRY3D_YLO] = wylo;
        width[BDRY3D_YHI] = wyhi;
        width[BDRY3D_ZLO] = wzlo;
        width[BDRY3D_ZHI] = wzhi;
    }

   
}
