#include "negative3d.h"
#include "utility.h"

namespace levelset {
        
    Negative3D::Negative3D(InitialFunc3D* g)
        : f(g), InitialFunc3D()
    {
        start[0] = 0; start[1] = 0;
        stop[0] = 1.; stop[1] = 1.;
        if (parameter) delete[] parameter;
    }

    double Negative3D::XYZ(const double s, const double t, const double u) const
    {
        return -f->XYZ(s,t,u);
    }

}
