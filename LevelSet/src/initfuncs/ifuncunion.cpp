#include "ifuncunion.h"
#include "utility.h"

namespace levelset {
        
    IFuncUnion::IFuncUnion(InitialFunc** f_in, const int fnum_in)
        : f(f_in), fnum(fnum_in), InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
    }

// NOTE: THIS ASSUMES THAT THE INTERIOR OF THE SHAPE HAS NEGATIVE LEVELSET VALUE
    double IFuncUnion::XY(const double s, const double t) const
    {
        double ans = f[0]->XY(s,t);
        for (int i=1; i<fnum; ++i)
            ans = min(ans,f[i]->XY(s,t));
        return ans;
    }

}
