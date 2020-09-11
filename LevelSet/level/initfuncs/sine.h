/*************************************************
    sine.h

    $Header: sine.h,v 1.1 2000/05/31 11:08:59 chopp Exp $

    $Log:       sine.h,v $
    * Revision 1.1  2000/05/31  11:08:59  11:08:59  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __SINE_H__
#define __SINE_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Sine : public InitialFunc 
    {
        enum {XNormal=0, YNormal, Hyperplane, Amp, Freq, Offset};
   
    public:

        Sine(const double xn = 0., const double yn = 1., const double b = 0.,
             const double a = 0., const double f = 1., const double off = 0.);
    Sine(InputParams *params) : InitialFunc() {SetParams(params);} 

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}
#endif
