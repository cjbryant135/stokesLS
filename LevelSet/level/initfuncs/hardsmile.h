/*************************************************
    hardsmile.h

    $Header: hardsmile.h,v 1.1 99/02/04 14:23:01 chopp Exp $

    $Log:       hardsmile.h,v $
    * Revision 1.1  99/02/04  14:23:01  14:23:01  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.3  98/03/04  14:59:57  14:59:57  chopp (David Chopp)
    * *** empty log message ***
    * 
    * Revision 1.3  98/03/02  12:58:31  12:58:31  chopp (David Chopp)
    * *** empty log message ***
    * 
    *************************************************/

#ifndef __HARDSMILE_H__
#define __HARDSMILE_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Hardsmile : public InitialFunc 
    {
        enum {Inner=0, Outer};
   
    public:

    Hardsmile(const double irad, const double orad) : InitialFunc(irad,orad) 
        {start = 0; stop = 2.*M_PI*orad;}
    Hardsmile(InputParams *params) : InitialFunc() {SetParams(params);} 

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}

#endif
