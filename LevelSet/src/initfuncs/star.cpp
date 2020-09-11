/*************************************************
    star.cpp

    $Header: star.cpp,v 2.1 99/01/06 13:59:53 chopp Exp $

    $Log:       star.cpp,v $
Revision 2.1  99/01/06  13:59:53  13:59:53  chopp (David Chopp)
*** none ***

Revision 1.3  98/03/04  14:59:33  14:59:33  chopp (David Chopp)
*** empty log message ***

Revision 1.3  98/03/02  12:57:38  12:57:38  chopp (David Chopp)
Updated to use InputParameters for initialization

Revision 1.2  97/12/18  11:27:45  11:27:45  chopp (David Chopp)
*** empty log message ***

Revision 1.1  97/12/04  10:06:29  10:06:29  chopp (David Chopp)
Initial revision

*************************************************/

#include "star.h"

namespace levelset {
        
    void Star::SetParams(InputParams *params)
    {
        if (parameter) delete[] parameter;
        parameter = new double[3];
        parameter[NPoints] = params->GetIntParam("star points");
        parameter[Inner] = params->GetDoubleParam("inner radius");
        parameter[Outer] = params->GetDoubleParam("outer radius");
    }

    double Star::X(const double t) const
    {
        double s = t;
        while (s < 0) s += stop;
        while (s > stop) s -= stop;
        return ((parameter[Inner]+parameter[Outer])/2.
                +(parameter[Outer]-parameter[Inner])
                *cos(parameter[NPoints]*s)/2.)*cos(s);
    }

    double Star::Y(const double t) const
    {
        double s = t;
        while (s < 0) s += stop;
        while (s > stop) s -= stop;
        return ((parameter[Inner]+parameter[Outer])/2.
                +(parameter[Outer]-parameter[Inner])
                *cos(parameter[NPoints]*s)/2.)*sin(s);
    }

/*
  double Star::XY(const double s, const double t) const
  {
  double r = sqrt(s*s+t*t);
  double theta = atan2(t,s);
  return 1.-r/((parameter[Inner]+parameter[Outer])/2.
  +(parameter[Outer]-parameter[Inner])
  *cos(parameter[NPoints]*theta)/2.);
  }
*/

}
