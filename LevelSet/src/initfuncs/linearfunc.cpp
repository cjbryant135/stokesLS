#include "linearfunc.h"
#include "utility.h"
#define A parameter[Normx]
#define B parameter[Normy]
#define C parameter[Zero]

namespace levelset {

    LinearFunc::LinearFunc(const double nx, const double ny, const double z)
        : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = nx;
        B = ny;
        C = z;
    }

    void LinearFunc::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = param->GetDoubleParam("LinearFunc -- X Normal");
        B = param->GetDoubleParam("LinearFunc -- Y Normal");
        C = param->GetDoubleParam("LinearFunc -- Constant");
    }

    double LinearFunc::X(const double t) const
    {
        return A;
    }

    double LinearFunc::Y(const double t) const
    {
        return A;
    }

    double LinearFunc::XY(const double s, const double t) const
    {
        return s*A+t*B+C;
    }

}
