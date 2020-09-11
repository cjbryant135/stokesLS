#include "negative.h"
#include "utility.h"

namespace levelset {
        
    Negative::Negative(InitialFunc* g)
        : f(g), InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
    }

    double Negative::XY(const double s, const double t) const
    {
        return -f->XY(s,t);
    }

}
