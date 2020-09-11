/*************************************************
    star.h

    $Header: star.h,v 1.1 99/02/04 14:23:07 chopp Exp $

    $Log:       star.h,v $
    * Revision 1.1  99/02/04  14:23:07  14:23:07  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.4  98/05/18  14:06:43  14:06:43  chopp (David Chopp)
    * *** empty log message ***
    * 
    * Revision 1.3  98/03/04  15:00:09  15:00:09  chopp (David Chopp)
    * *** empty log message ***
    * 
    * Revision 1.3  98/03/02  12:58:48  12:58:48  chopp (David Chopp)
    * Updated to use InputParameters for initialization
    * 
    * Revision 1.2  97/12/18  11:28:10  11:28:10  chopp (David Chopp)
    * *** empty log message ***
    * 
    * Revision 1.1  97/12/04  10:29:45  10:29:45  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __STAR_H__
#define __STAR_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Star : public InitialFunc 
    {
        enum {Inner=0, Outer, NPoints};
   
    public:

    Star(const int pts, const double irad, const double orad)
        : InitialFunc(irad,orad,pts) 
        {start = 0; stop = 2.*M_PI;}

    Star(InputParams *params) : InitialFunc(params) 
        {start = 0; stop = 2.*M_PI; SetParams(params);}

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
//   virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };

}
#endif
