/*************************************************
    hardsmile.cpp

    $Header: hardsmile.cpp,v 2.3 99/01/06 13:59:51 chopp Exp $

    $Log:       hardsmile.cpp,v $
Revision 2.3  99/01/06  13:59:51  13:59:51  chopp (David Chopp)
*** none ***

Revision 1.4  98/03/04  14:59:19  14:59:19  chopp (David Chopp)
*** empty log message ***

Revision 1.4  98/03/02  12:56:31  12:56:31  chopp (David Chopp)
Changed input so that problems could be switched without recompiling


*************************************************/

#include "hardsmile.h"
#include "utility.h"
#define A parameter[Inner]
#define B parameter[Outer]

namespace levelset {

    void Hardsmile::SetParams(InputParams *param)
    {
        if (parameter) delete[] parameter;
        parameter = new double[2];
        A = param->GetDoubleParam("Inner Radius");
        B = param->GetDoubleParam("Outer Radius");
    }

    double Hardsmile::X(const double t) const
    {
        double s = t;
        while (s < 0) s += stop;
        while (s > stop) s -= stop;
        if (s <= M_PI*parameter[Inner]) {
            return parameter[Inner]*cos(s/parameter[Inner]);
        } else {
            s -= M_PI*parameter[Inner];
            if (s <= M_PI*(parameter[Outer]-parameter[Inner])/2.) {
                return -(parameter[Inner]+parameter[Outer])/2
                    +(parameter[Outer]-parameter[Inner])/2
                    *cos(s*2/(parameter[Outer]-parameter[Inner]));
            } else {
                s -= M_PI*(parameter[Outer]-parameter[Inner])/2.;
                if (s <= M_PI*parameter[Outer]) {
                    return -parameter[Outer]*cos(s/parameter[Outer]);
                } else {
                    s -= M_PI*parameter[Outer];
                    return (parameter[Inner]+parameter[Outer])/2
                        +(parameter[Outer]-parameter[Inner])/2
                        *cos(s*2/(parameter[Outer]-parameter[Inner]));
                }
            }
        }
    }

    double Hardsmile::Y(const double t) const
    {
        double s = t;
        while (s < 0) s += stop;
        while (s > stop) s -= stop;
        if (s <= M_PI*parameter[Inner]) {
            return -parameter[Inner]*sin(s/parameter[Inner]);
        } else {
            s -= M_PI*parameter[Inner];
            if (s <= M_PI*(parameter[Outer]-parameter[Inner])/2.) {
                return (parameter[Outer]-parameter[Inner])/2
                    *sin(s*2/(parameter[Outer]-parameter[Inner]));
            } else {
                s -= M_PI*(parameter[Outer]-parameter[Inner])/2.;
                if (s <= M_PI*parameter[Outer]) {
                    return -parameter[Outer]*sin(s/parameter[Outer]);
                } else {
                    s -= M_PI*parameter[Outer];
                    return (parameter[Outer]-parameter[Inner])/2
                        *sin(s*2/(parameter[Outer]-parameter[Inner]));
                }
            }
        }
    }


    double Hardsmile::XY(const double s, const double t) const
    {
        double as = fabs(s);
        if (t > 0) {
            return (parameter[Outer]-parameter[Inner])/2
                -sqrt(sqr(as-(parameter[Inner]+parameter[Outer])/2)+t*t);
        } else {
            return (parameter[Outer]-parameter[Inner])/2
                -fabs(sqrt(as*as+t*t)-(parameter[Outer]+parameter[Inner])/2);
        }             
    }

}
