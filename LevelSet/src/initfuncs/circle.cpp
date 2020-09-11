/*************************************************
    circle.cpp

    $Header: /stash/chopp/cvsroot/extend/circle.cpp,v 1.1.1.1 2003/05/09 17:30:39 chopp Exp $

    $Log: circle.cpp,v $
    Revision 1.1.1.1  2003/05/09 17:30:39  chopp
    Base code, this code may or may not work, but it did at one time.


    Revision 1.1  2003/05/08 19:11:28  chopp
    *** empty log message ***

Revision 1.1  2000/06/01  11:58:51  11:58:51  chopp (David Chopp)
Initial revision

*************************************************/

#include "circle.h"
#include "utility.h"
#define A parameter[XCenter]
#define B parameter[YCenter]
#define C parameter[Radius]

namespace levelset {

    Circle::Circle(const  double xc, const double yc, const double radius) 
        : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = xc;
        B = yc;
        C = radius;
    }

    void Circle::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = param->GetDoubleParam("X-Center");
        B = param->GetDoubleParam("Y-Center");
        C = param->GetDoubleParam("Radius");
    }

    double Circle::X(const double t) const
    {
        return A+C*cos(t);
    }

    double Circle::Y(const double t) const
    {
        return B+C*sin(t);
    }

    double Circle::XY(const double s, const double t) const
    {
        return sqrt(sqr(s-A)+sqr(t-B))-C;
    }

}
