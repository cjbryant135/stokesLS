/*************************************************
    mirroredcubic.h
	 Author: Colton Bryant 
    Initially created 06/11/20 
	 *************************************************/

#ifndef __MIRROREDCUBIC_H__
#define __MIRROREDCUBIC_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class MirroredCubic : public InitialFunc 
    {
        enum {A=0, B, C, D};
   
    public:
		
	    MirroredCubic(const double a = 0., const double b = 0., const double c = 0., const double d = 0.);

   	 MirroredCubic(InputParams *params) : InitialFunc(params) 
		 	{start = -1.; stop = 1., SetParams(params);} //NOTE: Here I assumed xmin=-1 and xmax=1

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
//   virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}
#endif
