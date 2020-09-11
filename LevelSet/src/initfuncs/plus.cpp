/*************************************************
    plus.cpp

    $Header: plus.cpp,v 2.1 99/01/06 13:59:52 chopp Exp $

    $Log:       plus.cpp,v $
Revision 2.1  99/01/06  13:59:52  13:59:52  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:35  14:59:35  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:57:57  12:57:57  chopp (David Chopp)
Initial revision


*************************************************/

#include "plus.h"
#include "utility.h"
#define A parameter[Period]
#define B parameter[Width]
#define C parameter[Length]

namespace levelset {
        
    void Plus::SetParams(InputParams *params)
    {
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = params->GetDoubleParam("Period");
        B = params->GetDoubleParam("Width");
        C = params->GetDoubleParam("Length");
    }
 
/*
  double Plus::X(const double t) const
  {
  return ((parameter[Inner]+parameter[Outer])/2.
  +(parameter[Outer]-parameter[Inner])
  *cos(parameter[NPoints]*s)/2.)*cos(s);
  }

  double Plus::Y(const double t) const
  {
  double s = t;
  while (s < 0) s += stop;
  while (s > stop) s -= stop;
  return ((parameter[Inner]+parameter[Outer])/2.
  +(parameter[Outer]-parameter[Inner])
  *cos(parameter[NPoints]*s)/2.)*sin(s);
  }
*/
    double Plus::BaseFunc(const double s, const double t) const
    {
        double x = max(fabs(s),fabs(t));
        double y = min(fabs(s),fabs(t));
        if (x <= B/2.) {
            return sqrt(sqr(x-B/2.)+sqr(y-B/2.));
        }
        else if (x <= C+B/2.) {
            if (y <= B/2.) {
                return min(B/2.-y, C+B/2.-x);
            }
            else
                return B/2.-y;
        }
        else if (y <= B/2.) {
            return C+B/2.-x;
        }
        else
            return -sqrt(sqr(x-C-B/2.)+sqr(y-B/2.));
    }
   
    double Plus::XY(const double s, const double t) const
    {
        return max(BaseFunc(s-(A+C)/2.,t-(A-C)/2.),
                   BaseFunc(s-(A-C)/2.,t+(A+C)/2.),
                   BaseFunc(s+(A-C)/2.,t-(A+C)/2.),
                   BaseFunc(s+(A+C)/2.,t+(A-C)/2.));
    }

}



