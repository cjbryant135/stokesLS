/*************************************************
    wobble.h

    $Header: wobble.h,v 1.1 99/02/04 14:23:05 chopp Exp $

    $Log:       wobble.h,v $
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

#ifndef __WOBBLE_H__
#define __WOBBLE_H__

#include <math.h>
#include "initialfunc3d.h"

namespace levelset {
        
    class Wobble : public InitialFunc3D
    {
        enum {Radius=0, Length, Epsilon, Num_Of_Params};
   
    public:

        Wobble(const double rad, const double len, const double eps);
    Wobble(InputParams *params) : InitialFunc3D() {SetParams(params);} 

        virtual double XYZ(const double s, const double t, const double u) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}
#endif
