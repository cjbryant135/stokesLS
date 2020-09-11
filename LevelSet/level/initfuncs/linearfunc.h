#ifndef __LINEARFUNC_H__
#define __LINEARFUNC_H__

#include <math.h>
#include "../initialfunc.h"

namespace levelset {
        
    class LinearFunc : public InitialFunc 
    {
        enum {Normx=0, Normy, Zero};
   
    public:

        LinearFunc(const double nx = 1., const double ny = 0., const double z = 0.);
    LinearFunc(InputParams *params) : InitialFunc() {SetParams(params);} 

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}
#endif
