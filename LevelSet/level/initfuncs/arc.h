/*************************************************
 arc.h
 
 $Header: arc.h,v 1.1 99/02/04 14:23:03 chopp Exp $
 
 $Log:  arc.h,v $
 * Revision 1.1  99/02/04  14:23:03  14:23:03  chopp (David Chopp)
 * Initial revision
 * 
 *************************************************/

/*************************************************
 arc.h
 
 $Header: arc.h,v 1.1 99/02/04 14:23:03 chopp Exp $
 
 $Log:  arc.h,v $
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

#ifndef __ARC_H__
#define __ARC_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Arc : public InitialFunc 
    {
        enum {Height=0, Width};
                        
    public:
                        
    Arc(const int p, const double w, const double l)
        : InitialFunc(p,w,l) {}
    Arc(InputParams *params) : InitialFunc() {SetParams(params);}
                        
        virtual double XY(const double s, const double t) const;
                        
    protected:
                        
        virtual void SetParams(InputParams *params);
    };
        
}

#endif
