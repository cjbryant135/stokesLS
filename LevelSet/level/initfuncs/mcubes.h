/*************************************************
    mcubes.h

    $Header: mcubes.h,v 1.1 99/02/04 14:23:03 chopp Exp $

    $Log:       mcubes.h,v $
    * Revision 1.1  99/02/04  14:23:03  14:23:03  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  98/03/04  15:00:22  15:00:22  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  98/03/02  12:59:20  12:59:20  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __MCUBES_H__
#define __MCUBES_H__

#include <math.h>
#include "initialfunc3d.h"

namespace levelset {
        
    class Mcubes : public InitialFunc3D
    {
        int num;
   
    public:

    Mcubes(const int p, const double w, const double l)
        : InitialFunc3D(p,w,l) {}
    Mcubes(InputParams *params) : InitialFunc3D() {SetParams(params);}

        virtual double XYZ(const double s, const double t, const double u) const;

    protected:

        virtual void SetParams(InputParams *params);
    };

}

#endif
