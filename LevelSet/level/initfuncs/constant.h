#ifndef __CONSTANT_H__
#define __CONSTANT_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Constant : public InitialFunc 
    {
        enum {Value=0};
   
    public:

        Constant(const double val = 0.);
    Constant(InputParams *params) : InitialFunc() {SetParams(params);} 

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}
#endif
