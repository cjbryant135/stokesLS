/*************************************************
    wobble.cpp

    $Header: wobble.cpp,v 2.1 99/01/06 13:59:55 chopp Exp $

    $Log:       wobble.cpp,v $
Revision 2.1  99/01/06  13:59:55  13:59:55  chopp (David Chopp)
*** none ***

Revision 1.4  98/03/04  14:59:24  14:59:24  chopp (David Chopp)
*** empty log message ***

Revision 1.4  98/03/02  12:56:57  12:56:57  chopp (David Chopp)
Updated to use InputParameters for initialization


*************************************************/

#include "wobble.h"
#include "utility.h"
#define A parameter[Radius]
#define B parameter[Length]
#define C parameter[Epsilon]

namespace levelset {
        
    Wobble::Wobble(const double rad, const double len, const double eps)
        : InitialFunc3D()
    {
        if (parameter) delete[] parameter;
        parameter = new double[Num_Of_Params];
        A = rad;
        B = len;
        C = eps;
    }

    void Wobble::SetParams(InputParams *param)
    {
        if (parameter) delete[] parameter;
        parameter = new double[Num_Of_Params];
        A = param->GetDoubleParam("Radius",1.);
        B = param->GetDoubleParam("Period Length",1.);
        C = param->GetDoubleParam("Epsilon",0.01,"Size of perturbation");
    }

    double Wobble::XYZ(const double s, const double t, const double u) const
    {
        return A-C*cos(s*2*M_PI/B)-sqrt(t*t+u*u);
    }

}
