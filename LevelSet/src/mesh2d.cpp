#include "mesh2d.h"
#include "utility.h"

namespace levelset {
        
    char Mesh2D::IsCrossing(const double value, const double y0, const double y1,
                            double& f)
    {
        if (value-y0 >= 0 && value-y0 <= y1-y0) {
            f = frac(value, y0, y1);
            return true;
        } else
            return false;
    }

}
