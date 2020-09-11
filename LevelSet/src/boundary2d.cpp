/*************************************************
    boundary2d.cpp

    $Header: boundary2d.cpp,v 2.3 99/01/06 13:59:03 chopp Exp $

    $Log:       boundary2d.cpp,v $
Revision 2.3  99/01/06  13:59:03  13:59:03  chopp (David Chopp)
*** none ***

Revision 1.1  97/12/04  10:05:20  10:05:20  chopp (David Chopp)
Initial revision

*************************************************/

#include "boundary2d.h"

namespace levelset {

    Boundary2D::Boundary2D(const int wxlo, const int wxhi,
                           const int wylo, const int wyhi)
    {
        width[BDRY2D_XLO] = wxlo;
        width[BDRY2D_XHI] = wxhi;
        width[BDRY2D_YLO] = wylo;
        width[BDRY2D_YHI] = wyhi;
    }

   
}
