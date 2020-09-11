/*************************************************
    ovoid.cpp

    $Header: ovoid.cpp,v 2.1 99/01/06 13:59:55 chopp Exp $

    $Log:       ovoid.cpp,v $
Revision 2.1  99/01/06  13:59:55  13:59:55  chopp (David Chopp)
*** none ***

Revision 1.4  98/03/04  14:59:24  14:59:24  chopp (David Chopp)
*** empty log message ***

Revision 1.4  98/03/02  12:56:57  12:56:57  chopp (David Chopp)
Updated to use InputParameters for initialization


*************************************************/

#include "ovoid.h"
#include "utility.h"
#define A parameter[XCenter]
#define B parameter[YCenter]
#define C parameter[ZCenter]
#define D parameter[XWidth]
#define E parameter[YWidth]
#define F parameter[ZWidth]
#define G parameter[Angle1]
#define H parameter[Angle2]

namespace levelset {
        
    Ovoid::Ovoid(const  double xc, const double yc, const double zc,
                 const double xw, const double yw, const double zw,
                 const double angle1, const double angle2)
        : InitialFunc3D()
    {
        if (parameter) delete[] parameter;
        parameter = new double[Num_Of_Params];
        A = xc;
        B = yc;
        C = zc;
        D = xw;
        E = yw;
        F = zw;
        G = angle1;
        H = angle2;
    }

    void Ovoid::SetParams(InputParams *param)
    {
        if (parameter) delete[] parameter;
        parameter = new double[Num_Of_Params];
        A = param->GetDoubleParam("X-Center",0.,"X Center");
        B = param->GetDoubleParam("Y-Center",0.,"Y Center");
        C = param->GetDoubleParam("Z-Center",0.,"Z Center");
        D = param->GetDoubleParam("X Width",1.,"X Diameter");
        E = param->GetDoubleParam("Y Width",1.,"Y Diameter");
        F = param->GetDoubleParam("Z Width",1.,"Z Diameter");
        G = param->GetDoubleParam("Z Rotation",0.,"Tilt Angle");
        H = param->GetDoubleParam("X-Y Rotation",0.,"Plane Angle");
    }

    double Ovoid::XYZ(const double s, const double t, const double u) const
    {
        double as = (s-A)*cos(G)+sin(G)*((t-B)*cos(H)+(u-C)*sin(H));
        double at = -(s-A)*sin(G)+cos(G)*((t-B)*cos(H)+(u-C)*sin(H));
        double au = -(t-B)*sin(H)+(u-C)*cos(H);
        return 1-sqrt(sqr(as/D)+sqr(at/E)+sqr(au/F));
    }

}
