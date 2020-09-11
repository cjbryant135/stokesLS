/*************************************************
    eltype.h

    $Header: eltype.h,v 1.1 99/09/20 11:36:07 chopp Exp $

    $Log:       eltype.h,v $
    * Revision 1.1  99/09/20  11:36:07  11:36:07  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  99/02/26  14:26:34  14:26:34  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __ELTYPE_H__
#define __ELTYPE_H__

#include <math.h>
#include <float.h>

namespace blockmatrix {

    inline double IsZero(const double a)
    {
        return fabs(a) < 10.0 * DBL_MIN;
    }

    inline double IsOne(const double a)
    {
        return fabs(a - 1.0) < 10.0 * DBL_EPSILON;
    }

    struct Element 
    {
    public:

        size_t Pos;
        double Val;
    };
    typedef struct Element ElType;

    typedef int myBool;

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

    const ElType ZeroEl = {0, 0.};

}

#endif
