/*************************************************
    sheet.cpp

    $Header: sheet.cpp,v 2.1 99/01/06 13:59:51 chopp Exp $

    $Log:       sheet.cpp,v $
Revision 2.1  99/01/06  13:59:51  13:59:51  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:49  14:59:49  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:58:14  12:58:14  chopp (David Chopp)
Initial revision


*************************************************/

#include "sheet.h"
#include "utility.h"
#include <limits.h>
#define A parameter[Amplitude]
#define B parameter[XWidth]
#define C parameter[YWidth]
#define D parameter[XCycles]
#define E parameter[YCycles]

namespace levelset {
        
    void Sheet::SetParams(InputParams *params)
    {
        if (parameter) delete[] parameter;
        parameter = new double[5];
        A = params->GetDoubleParam("Amplitude");
        B = params->GetDoubleParam("xmax");
        C = params->GetDoubleParam("ymax");
        D = params->GetDoubleParam("X-cycles");
        E = params->GetDoubleParam("Y-cycles");
    }

    double Sheet::XYZ(const double s, const double t, const double u) const
    {
        return u-A*cos(D*s*M_PI/B)*cos(E*t*M_PI/C);
    }



}





