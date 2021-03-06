/*************************************************
    cup.h

    $Header: cup.h,v 1.1 99/02/04 14:22:56 chopp Exp $

    $Log:       cup.h,v $
    * Revision 1.1  99/02/04  14:22:56  14:22:56  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  98/03/04  15:00:18  15:00:18  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  98/03/02  12:59:13  12:59:13  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __CUP_H__
#define __CUP_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Cup : public InitialFunc 
    {
        enum {DotWidth=0, InsideC, OutsideC};
   
    public:

    Cup(const double p, const double w, const double l)
        : InitialFunc(p,w,l) {} 
    Cup(InputParams *params) : InitialFunc() {SetParams(params);}

        double CFunc(const double s, const double t) const;
        double DotFunc(const double s, const double t) const;
   
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}

#endif
