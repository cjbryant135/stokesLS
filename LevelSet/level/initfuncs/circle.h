/*************************************************
    circle.h

    $Header: /stash/chopp/cvsroot/extend/circle.h,v 1.1.1.1 2003/05/09 17:30:39 chopp Exp $

    $Log: circle.h,v $
    Revision 1.1.1.1  2003/05/09 17:30:39  chopp
    Base code, this code may or may not work, but it did at one time.


    Revision 1.1  2003/05/08 19:15:16  chopp
    Base code, may or may not work

    * Revision 1.1  2000/06/01  12:00:02  12:00:02  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __CIRCLE_H__
#define __CIRCLE_H__

#include <math.h>
#include "../initialfunc.h"

namespace levelset {

    class Circle : public InitialFunc 
    {
        enum {XCenter=0, YCenter, Radius};
   
    public:

    Circle(void) : InitialFunc() {}
        Circle(const double xc, const double yc, const double radius);
    Circle(InputParams *params) : InitialFunc() {SetParams(params);} 

        virtual double X(const double t) const;
        virtual double Y(const double t) const;
        virtual double XY(const double s, const double t) const;

    protected:
        virtual void SetParams(InputParams *params);
    };
        
}

#endif
