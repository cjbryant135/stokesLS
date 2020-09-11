/*************************************************
    line.cpp

    $Header: line.cpp,v 1.1 2000/05/31 11:07:56 chopp Exp $

    $Log:       line.cpp,v $
Revision 1.1  2000/05/31  11:07:56  11:07:56  chopp (David Chopp)
Initial revision

*************************************************/

#include "line.h"
#include "utility.h"
#define NX parameter[XNormal]
#define NY parameter[YNormal]
#define B parameter[Hyperplane]

namespace levelset {

    Line::Line(const  double xn, const double yn, const double b)
        : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        NX = xn;
        NY = yn;
        B = b;
    }

    void Line::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[5];
        NX = param->GetDoubleParam("X-Normal",0.,"X-component of line normal");
        NY = param->GetDoubleParam("Y-Normal",1.,"Y-component of line normal");
        B = param->GetDoubleParam("Hyperplane",0., 
                                  "Line formed by (X-normal,Y-normal).(x,y)=Hyperplane");
    }

    double Line::X(const double t) const
    {
        return t;
    }

    double Line::Y(const double t) const
    {
        return 0.;
    }

    double Line::XY(const double s, const double t) const
    {
        return (NX*s+NY*t)/sqrt(NX*NX+NY*NY)-B;
    }

}
