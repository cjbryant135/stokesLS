#ifndef __NEGATIVE3D_H__
#define __NEGATIVE3D_H__

#include <math.h>
#include "initialfunc3d.h"

namespace levelset {
        
    class Negative3D : public InitialFunc3D
    {
        InitialFunc3D* f;
        
    public:

        Negative3D(InitialFunc3D* g);

        virtual double XYZ(const double s, const double t, const double u) const;
    };

}
#endif
