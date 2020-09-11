/*************************************************
    oval.cpp

    $Header: /stash/chopp/cvsroot/extend/oval.cpp,v 1.1.1.1 2003/05/09 17:30:39 chopp Exp $

    $Log: oval.cpp,v $
    Revision 1.1.1.1  2003/05/09 17:30:39  chopp
    Base code, this code may or may not work, but it did at one time.


    Revision 1.1  2003/05/08 19:15:16  chopp
    Base code, may or may not work

Revision 1.1  2000/06/01  11:59:20  11:59:20  chopp (David Chopp)
Initial revision

*************************************************/

#include "oval.h"
#include "utility.h"
#define A parameter[XCenter]
#define B parameter[YCenter]
#define C parameter[XWidth]
#define D parameter[YWidth]
#define E parameter[Angle]

namespace levelset {
        
    Oval::Oval(const  double xc, const double yc, const double xw, const double yw,
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

    void Oval::SetParams(InputParams *param)
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

    double Oval::X(const double t) const
    {
        return parameter[XCenter] + parameter[XWidth]*parameter[YWidth]
            *cos(t+parameter[Angle])/sqrt(
                sqr(parameter[XWidth]*sin(t+parameter[Angle]))
                +sqr(parameter[YWidth]*cos(t+parameter[Angle])));
    }

    double Oval::Y(const double t) const
    {
        return parameter[YCenter] + parameter[XWidth]*parameter[YWidth]
            *sin(t+parameter[Angle])/sqrt(
                sqr(parameter[XWidth]*sin(t+parameter[Angle]))
                +sqr(parameter[YWidth]*cos(t+parameter[Angle])));
    }

    double Oval::XY(const double s, const double t) const
    {
        double as = (s-A)*cos(parameter[Angle])+(t-B)*sin(parameter[Angle]);
        double at = -(s-A)*sin(parameter[Angle])+(t-B)*cos(parameter[Angle]);
        return 1-sqrt(sqr(as/parameter[XWidth])+sqr(at/parameter[YWidth]));
    }

}
