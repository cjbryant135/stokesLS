/*************************************************
    plus.h

    $Header: plus.h,v 1.1 99/02/04 14:23:06 chopp Exp $

    $Log:       plus.h,v $
    * Revision 1.1  99/02/04  14:23:06  14:23:06  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  98/03/04  15:00:15  15:00:15  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  98/03/02  12:59:05  12:59:05  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __PLUS_H__
#define __PLUS_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Plus : public InitialFunc 
    {
        enum {Period=0, Width, Length};
   
    public:

    Plus(const double p, const double w, const double l)
        : InitialFunc(p,w,l) {} 
    Plus(InputParams *params) : InitialFunc(params) {}

        double BaseFunc(const double s, const double t) const;
   
//   virtual double X(const double t) const;
//   virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}
#endif
