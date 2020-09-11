#include "heapelt3d.h"

namespace levelset {
        
    HeapElement3D::HeapElement3D(const int i, const int j, const int k,
                                 const double val) 
        : value(val)
    { ind[0] = i; ind[1] = j; ind[2] = k;}



}





