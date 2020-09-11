#include "wavyfunc.h"
#include "utility.h"
#define A parameter[Wavenum]
#define B parameter[Amplitude]
#define C parameter[Low]

namespace levelset {
        
    WavyFunc::WavyFunc(const double wv, const double a, const double lo)
        : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = wv;
        B = a;
        C = lo;
    }

    void WavyFunc::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[3];
        A = param->GetDoubleParam("WavyFunc -- Wave Number");
        B = param->GetDoubleParam("WavyFunc -- Amplitude");
        C = param->GetDoubleParam("WavyFunc -- Minimum");
    }

    double WavyFunc::X(const double t) const
    {
        return A;
    }

    double WavyFunc::Y(const double t) const
    {
        return A;
    }

    double WavyFunc::XY(const double s, const double t) const
    {
        double theta = atan2(t,s);
        return C+B*(cos(A*theta)+1);
    }

}
