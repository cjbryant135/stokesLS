/*************************************************
    diamond.h

    $Header: /stash/chopp/cvsroot/extend/diamond.h,v 1.1.1.1 2003/05/09 17:30:39 chopp Exp $

    $Log: diamond.h,v $
    Revision 1.1.1.1  2003/05/09 17:30:39  chopp
    Base code, this code may or may not work, but it did at one time.


    Revision 1.1  2003/05/08 19:15:16  chopp
    Base code, may or may not work

    * Revision 1.1  2000/06/22  11:28:13  11:28:13  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  2000/06/01  12:00:02  12:00:02  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __DIAMOND_H__
#define __DIAMOND_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {

    class Diamond : public InitialFunc 
    {
        enum {XCenter=0, YCenter, Radius};
   
    public:

        Diamond(const double xc, const double yc, const double radius);
    Diamond(InputParams *params) : InitialFunc() {SetParams(params);} 

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };
        
}

#endif
