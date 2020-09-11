/*************************************************
    dblrect.cpp

    $Header: dblrect.cpp,v 1.4 98/03/04 14:59:24 chopp Exp $

    $Log:       dblrect.cpp,v $
Revision 1.4  98/03/04  14:59:24  14:59:24  chopp (David Chopp)
*** empty log message ***

Revision 1.4  98/03/02  12:56:57  12:56:57  chopp (David Chopp)
Updated to use InputParameters for initialization


*************************************************/

#include "dblrect.h"
#include "utility.h"
#define A parameter[XCenter]
#define B parameter[YCenter]
#define C parameter[XWidth]
#define D parameter[YWidth]
#define E parameter[Angle]

namespace levelset {

    DblRect::DblRect(const  double xc, const double yc, const double xw, const double yw,
                     const double angle)
        : InitialFunc()
    {
        start = 0;
        stop = 4.;
        if (parameter) delete[] parameter;
        parameter = new double[5];
        A = xc;
        B = yc;
        C = xw;
        D = yw;
        E = angle;
    }

    void DblRect::SetParams(InputParams *param)
    {
        start = 0;
        stop = 4.;
        if (parameter) delete[] parameter;
        parameter = new double[5];
        A = param->GetDoubleParam("X-Center");
        B = param->GetDoubleParam("Y-Center");
        C = param->GetDoubleParam("X Width");
        D = param->GetDoubleParam("Y Width");
        E = param->GetDoubleParam("Rotation");
    }

    double DblRect::X(const double t) const
    {
        return parameter[XCenter] + parameter[XWidth]*parameter[YWidth]
            *cos(t+parameter[Angle])/sqrt(
                sqr(parameter[XWidth]*sin(t+parameter[Angle]))
                +sqr(parameter[YWidth]*cos(t+parameter[Angle])));
    }

    double DblRect::Y(const double t) const
    {
        return parameter[YCenter] + parameter[XWidth]*parameter[YWidth]
            *sin(t+parameter[Angle])/sqrt(
                sqr(parameter[XWidth]*sin(t+parameter[Angle]))
                +sqr(parameter[YWidth]*cos(t+parameter[Angle])));
    }

    double DblRect::XY(const double s, const double t) const
    {
        return min(parameter[YCenter]+parameter[YWidth]-fabs(t),
                   fabs(t)-parameter[YCenter]+parameter[YWidth],
                   parameter[XCenter]+parameter[XWidth]-fabs(s),
                   fabs(s)-parameter[XCenter]+parameter[XWidth]);
    }

}

