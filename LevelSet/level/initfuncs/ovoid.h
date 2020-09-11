/*************************************************
    ovoid.h

    $Header: ovoid.h,v 1.1 99/02/04 14:23:05 chopp Exp $

    $Log:       ovoid.h,v $
    * Revision 1.1  99/02/04  14:23:05  14:23:05  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.4  98/03/04  15:00:04  15:00:04  chopp (David Chopp)
    * *** empty log message ***
    * 
    * Revision 1.4  98/03/02  12:58:42  12:58:42  chopp (David Chopp)
    * Updated to use InputParameters for initialization
    * 
    *************************************************/

#ifndef __OVOID_H__
#define __OVOID_H__

#include <math.h>
#include "initialfunc3d.h"

namespace levelset {
        
    class Ovoid : public InitialFunc3D
    {
        enum {XCenter=0, YCenter, ZCenter, XWidth, YWidth, ZWidth, Angle1, Angle2,
              Num_Of_Params};
   
    public:

        Ovoid(const double xc, const double yc, const double zc,
              const double xw, const double yw, const double zw,
              const double angle1 = 0, const double angle2 = 0);
    Ovoid(InputParams *params) : InitialFunc3D() {SetParams(params);} 

        virtual double XYZ(const double s, const double t, const double u) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}
#endif
