#include "heavifunc.h"
#include "utility.h"
#define A parameter[Normx]
#define B parameter[Normy]
#define C parameter[Beta]
#define D parameter[Value1]
#define E parameter[Value2]

namespace levelset {
        
    HeavisideFunc::HeavisideFunc(const double nx, const double ny, const double b, const double v1,
                                 const double v2) : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[5];
        A = nx;
        B = ny;
        C = b;
        D = v1;
        E = v2;
    }

    void HeavisideFunc::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[5];
        A = param->GetDoubleParam("HeavisideFunc -- X Normal");
        B = param->GetDoubleParam("HeavisideFunc -- Y Normal");
        C = param->GetDoubleParam("HeavisideFunc -- Beta");
        D = param->GetDoubleParam("HeavisideFunc -- Low Value");
        E = param->GetDoubleParam("HeavisideFunc -- High Value");
    }

    double HeavisideFunc::X(const double t) const
    {
        return A;
    }

    double HeavisideFunc::Y(const double t) const
    {
        return A;
    }

    double HeavisideFunc::XY(const double s, const double t) const
    {
        return (A*s+B*t < C) ? D : E;
    }

}
