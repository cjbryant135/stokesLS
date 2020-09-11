/*************************************************
    plus2.cpp

    $Header: plus2.cpp,v 2.1 99/01/06 13:59:53 chopp Exp $

    $Log:       plus2.cpp,v $
Revision 2.1  99/01/06  13:59:53  13:59:53  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:44  14:59:44  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:58:03  12:58:03  chopp (David Chopp)
Initial revision


*************************************************/

#include "plus2.h"
#include "utility.h"
#define A parameter[Period]
#define B parameter[Width]
#define C parameter[Length]

namespace levelset {
        
    void Plus2::SetParams(InputParams *param)
    {
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = param->GetDoubleParam("Period");
        B = param->GetDoubleParam("Width");
        C = param->GetDoubleParam("Length");
    }


/*
  double Plus2::X(const double t) const
  {
  return ((parameter[Inner]+parameter[Outer])/2.
  +(parameter[Outer]-parameter[Inner])
  *cos(parameter[NPoints]*s)/2.)*cos(s);
  }

  double Plus2::Y(const double t) const
  {
  double s = t;
  while (s < 0) s += stop;
  while (s > stop) s -= stop;
  return ((parameter[Inner]+parameter[Outer])/2.
  +(parameter[Outer]-parameter[Inner])
  *cos(parameter[NPoints]*s)/2.)*sin(s);
  }
*/
    double Plus2::BaseFunc(const double s, const double t) const
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
   
    double Plus2::XY(const double s, const double t) const
    {
        return max(BaseFunc(s,t),C+B/2.-A-sqrt(sqr(fabs(s)-A)+sqr(fabs(t)-A)));
    }

}



