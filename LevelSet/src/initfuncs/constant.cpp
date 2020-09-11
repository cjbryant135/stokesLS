#include "constant.h"
#include "utility.h"
#define A parameter[Value]

namespace levelset {

    Constant::Constant(const double val)
        : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[2];
        A = val;
    }

    void Constant::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[2];
        A = param->GetDoubleParam("Constant Value",0.,"Value of a constant");
    }

    double Constant::X(const double t) const
    {
        return A;
    }

    double Constant::Y(const double t) const
    {
        return A;
    }

    double Constant::XY(const double s, const double t) const
    {
        return A;
    }

}
