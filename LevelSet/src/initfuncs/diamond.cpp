/*************************************************
    diamond.cpp

    $Header: /stash/chopp/cvsroot/extend/diamond.cpp,v 1.1.1.1 2003/05/09 17:30:39 chopp Exp $

    $Log: diamond.cpp,v $
    Revision 1.1.1.1  2003/05/09 17:30:39  chopp
    Base code, this code may or may not work, but it did at one time.


    Revision 1.1  2003/05/08 19:15:16  chopp
    Base code, may or may not work

Revision 1.1  2000/06/22  11:28:05  11:28:05  chopp (David Chopp)
Initial revision

*************************************************/

/*************************************************
    diamond.cpp

    $Header: /stash/chopp/cvsroot/extend/diamond.cpp,v 1.1.1.1 2003/05/09 17:30:39 chopp Exp $

    $Log: diamond.cpp,v $
    Revision 1.1.1.1  2003/05/09 17:30:39  chopp
    Base code, this code may or may not work, but it did at one time.


    Revision 1.1  2003/05/08 19:15:16  chopp
    Base code, may or may not work

Revision 1.1  2000/06/22  11:28:05  11:28:05  chopp (David Chopp)
Initial revision

Revision 1.1  2000/06/01  11:58:51  11:58:51  chopp (David Chopp)
Initial revision

*************************************************/

#include "diamond.h"
#include "utility.h"
#define A parameter[XCenter]
#define B parameter[YCenter]
#define C parameter[Radius]

namespace levelset {

    Diamond::Diamond(const  double xc, const double yc, const double radius) 
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

    void Diamond::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = param->GetDoubleParam("X-Center");
        B = param->GetDoubleParam("Y-Center");
        C = param->GetDoubleParam("Radius");
    }

    double Diamond::X(const double t) const
    {
        return A+C*cos(t);
    }

    double Diamond::Y(const double t) const
    {
        return B+C*sin(t);
    }

    double Diamond::XY(const double s, const double t) const
    {
        if (C-fabs(s-A)-fabs(t-B) > 0.) {
            return (C-fabs(s-A)-fabs(t-B))/sqrt(2.);
        }
        else if (fabs(fabs(s-A)-fabs(t-B)) < C) {
            return -fabs(C-fabs(s-A)-fabs(t-B))/sqrt(2.);
        }
        else return max(-sqrt(sqr(fabs(s-A)-C)+(t-B)*(t-B)),
                        -sqrt((s-A)*(s-A)+sqr(fabs(t-B)-C)));
    }

}
