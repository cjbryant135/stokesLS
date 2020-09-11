/*************************************************
    dbloval.h

    $Header: dbloval.h,v 1.1 99/02/04 14:22:58 chopp Exp $

    $Log:       dbloval.h,v $
    * Revision 1.1  99/02/04  14:22:58  14:22:58  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.4  98/03/04  15:00:04  15:00:04  chopp (David Chopp)
    * *** empty log message ***
    * 
    * Revision 1.4  98/03/02  12:58:42  12:58:42  chopp (David Chopp)
    * Updated to use InputParameters for initialization
    * 
    *************************************************/

#ifndef __DBLOVAL_H__
#define __DBLOVAL_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class DblOval : public InitialFunc 
    {
        enum {XCenter=0, YCenter, XWidth, YWidth, Angle};
   
    public:

        DblOval(const double xc, const double yc, const double xw, const double yw,
                const double angle = 0);
    DblOval(InputParams *params) : InitialFunc() {SetParams(params);} 

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}

#endif
