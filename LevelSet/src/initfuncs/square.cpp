/*************************************************
    square.cpp

    $Header: /stash/chopp/cvsroot/extend/square.cpp,v 1.1.1.1 2003/05/09 17:30:39 chopp Exp $

    $Log: square.cpp,v $
    Revision 1.1.1.1  2003/05/09 17:30:39  chopp
    Base code, this code may or may not work, but it did at one time.


    Revision 1.1  2003/05/08 19:15:16  chopp
    Base code, may or may not work

Revision 1.2  2000/06/09  22:05:57  22:05:57  chopp (David Chopp)
*** none ***

Revision 1.1  2000/06/01  11:59:25  11:59:25  chopp (David Chopp)
Initial revision

*************************************************/

#include "square.h"
#include "utility.h"
#define A parameter[XCenter]
#define B parameter[YCenter]
#define C parameter[Radius]

namespace levelset {
        
    Square::Square(const  double xc, const double yc, const double radius) 
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

    void Square::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = param->GetDoubleParam("X-Center");
        B = param->GetDoubleParam("Y-Center");
        C = param->GetDoubleParam("Radius");
    }

    double Square::X(const double t) const
    {
        return A+C*cos(t);
    }

    double Square::Y(const double t) const
    {
        return B+C*sin(t);
    }

    double Square::XY(const double s, const double t) const
    {
        double ans;
        if (fabs(s-A) < C) {
            if (fabs(t-B) < C) {
                ans = C-max(fabs(s-A),fabs(t-B));
            }
            else {
                ans = C-fabs(t-B);
            }
        }
        else {
            if (fabs(t-B) < C) {
                ans = C-fabs(s-A);
            }
            else {
                if (t-B >= C) {
                    if (s-A >= C) {
                        ans = -sqrt(sqr(t-B-C)+sqr(s-A-C));
                    }
                    else {
                        ans = -sqrt(sqr(t-B-C)+sqr(s-A+C));
                    }
                }
                else {
                    if (s-A >= C) {
                        ans = -sqrt(sqr(t-B+C)+sqr(s-A-C));
                    }
                    else {
                        ans = -sqrt(sqr(t-B+C)+sqr(s-A+C));
                    }
                }
            }
        }
        return ans;
    }
}

