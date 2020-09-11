/*************************************************
    snake.h

    $Header: snake.h,v 1.1 99/02/04 14:23:03 chopp Exp $

    $Log:       snake.h,v $
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

#ifndef __SNAKE_H__
#define __SNAKE_H__

#include <math.h>
#include "initialfunc3d.h"

namespace levelset {
        
    class Snake : public InitialFunc3D
    {
        enum{Length=0, Width, Direction};
        const int Num;
        double *snake[3];
    public:

    Snake(const int p, const double w, const double l)
        : InitialFunc3D(p,w,l), Num(1000) {}
    Snake(InputParams *params) : InitialFunc3D(), Num(1000) {SetParams(params);}
        ~Snake(void);

        virtual double XYZ(const double s, const double t, const double u) const;

    protected:

        virtual void SetParams(InputParams *params);
        virtual double XPath(const double s);
        virtual double YPath(const double s);
        virtual double ZPath(const double s);
    };

}
#endif
