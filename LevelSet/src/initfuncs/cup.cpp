/*************************************************
    cup.cpp

    $Header: cup.cpp,v 2.4 99/01/06 13:59:50 chopp Exp $

    $Log:       cup.cpp,v $
Revision 2.4  99/01/06  13:59:50  13:59:50  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:46  14:59:46  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:58:06  12:58:06  chopp (David Chopp)
Initial revision


*************************************************/

#include "cup.h"
#include "utility.h"
#define A parameter[DotWidth]
#define B parameter[InsideC]
#define C parameter[OutsideC]

namespace levelset {

    void Cup::SetParams(InputParams *param)
    {
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = param->GetDoubleParam("Dot Width");
        B = param->GetDoubleParam("C Inner Width");
        C = param->GetDoubleParam("C Outer Width");
    }


    double Cup::CFunc(const double s, const double t) const
    {
        double tt = fabs(t);
        if (s < -C) {
            if (tt > C) {
                return -sqrt(sqr(s+C)+sqr(tt-C));
            }
            else
                return s+C;
        }
        else {
            if (s < -B) {
                if (tt <= B) {
                    return min(-(B+s),s+C);
                }
                else {
                    if (tt <= C) {
                        return min(C-tt,s+C,sqrt(sqr(B+s)+sqr(tt-B)));
                    }
                    else
                        return C-tt;
                }
            }
            else {
                if (s < C) {
                    if (tt < B) {
                        return -min(B-tt,s+B);
                    }
                    else {
                        if (tt < C) {
                            return min(tt-B,C-tt,C-s);
                        }
                        else
                            return C-tt;
                    }
                }
                else {
                    if (tt < B) {
                        return -sqrt(sqr(s-C)+sqr(B-tt));
                    }
                    else {
                        if (tt < C) {
                            return C-s;
                        }
                        else
                            return -sqrt(sqr(s-C)+sqr(tt-C));
                    }
                }
            }
        }
    }

    double Cup::DotFunc(const double s, const double t) const
    {
        double ss = fabs(s);
        double tt = fabs(t);
        if (ss < A) {
            if (tt < A) {
                return min(A-ss,A-tt);
            }
            else
                return A-tt;
        }
        else {
            if (tt < A) {
                return A-ss;
            }
            else
                return -sqrt(sqr(A-ss)+sqr(A-tt));
        }
    }


    double Cup::XY(const double s, const double t) const
    {
        return max(DotFunc(s,t),CFunc(s,t));
    }


}






