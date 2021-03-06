/*************************************************
    plus2.h

    $Header: plus2.h,v 1.1 99/02/04 14:23:07 chopp Exp $

    $Log:       plus2.h,v $
    * Revision 1.1  99/02/04  14:23:07  14:23:07  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  98/03/04  15:00:17  15:00:17  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  98/03/02  12:59:10  12:59:10  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __PLUS2_H__
#define __PLUS2_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Plus2 : public InitialFunc 
    {
        enum {Period=0, Width, Length};
   
    public:

    Plus2(const double p, const double w, const double l)
        : InitialFunc(p,w,l) {} 
    Plus2(InputParams *params) : InitialFunc(params) {}

        double BaseFunc(const double s, const double t) const;
   
//   virtual double X(const double t) const;
//   virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}

#endif
