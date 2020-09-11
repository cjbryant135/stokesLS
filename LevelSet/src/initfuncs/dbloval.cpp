/*************************************************
    dbloval.cpp

    $Header: dbloval.cpp,v 1.4 98/03/04 14:59:24 chopp Exp $

    $Log:       dbloval.cpp,v $
Revision 1.4  98/03/04  14:59:24  14:59:24  chopp (David Chopp)
*** empty log message ***

Revision 1.4  98/03/02  12:56:57  12:56:57  chopp (David Chopp)
Updated to use InputParameters for initialization


*************************************************/

#include "dbloval.h"
#include "utility.h"
#define A parameter[XCenter]
#define B parameter[YCenter]
#define C parameter[XWidth]
#define D parameter[YWidth]
#define E parameter[Angle]

namespace levelset {
        
    DblOval::DblOval(const  double xc, const double yc, const double xw, const double yw,
                     const double angle)
        : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[5];
        A = xc;
        B = yc;
        C = xw;
        D = yw;
        E = angle;
    }

    void DblOval::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[5];
        A = param->GetDoubleParam("X-Center");
        B = param->GetDoubleParam("Y-Center");
        C = param->GetDoubleParam("X Width");
        D = param->GetDoubleParam("Y Width");
        E = param->GetDoubleParam("Rotation");
    }

    double DblOval::X(const double t) const
    {
        return parameter[XCenter] + parameter[XWidth]*parameter[YWidth]
            *cos(t+parameter[Angle])/sqrt(
                sqr(parameter[XWidth]*sin(t+parameter[Angle]))
                +sqr(parameter[YWidth]*cos(t+parameter[Angle])));
    }

    double DblOval::Y(const double t) const
    {
        return parameter[YCenter] + parameter[XWidth]*parameter[YWidth]
            *sin(t+parameter[Angle])/sqrt(
                sqr(parameter[XWidth]*sin(t+parameter[Angle]))
                +sqr(parameter[YWidth]*cos(t+parameter[Angle])));
    }

    double DblOval::XY(const double s, const double t) const
    {
        double as = s*cos(parameter[Angle])+t*sin(parameter[Angle]);
        double at = -s*sin(parameter[Angle])+t*cos(parameter[Angle]);
        return max(1-sqrt(sqr(as/parameter[XWidth])+sqr((at-parameter[YCenter])
                                                        /parameter[YWidth])),
                   1-sqrt(sqr(as/parameter[XWidth])+sqr((at+parameter[YCenter])
                                                        /parameter[YWidth])));
    }

}
