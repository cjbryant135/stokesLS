/*************************************************
    mdisks.h

    $Header: mdisks.h,v 1.1 99/02/04 14:23:03 chopp Exp $

    $Log:       mdisks.h,v $
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

#ifndef __MDISKS_H__
#define __MDISKS_H__

#include <math.h>
#include "initialfunc.h"

namespace levelset {
        
    class Mdisks : public InitialFunc 
    {
        int num;
   
    public:

    Mdisks(void) : InitialFunc() {}
    Mdisks(const int p, const double w, const double l)
        : InitialFunc(p,w,l) {}
    Mdisks(InputParams *params) : InitialFunc() {SetParams(params);}
    Mdisks(const int n, const double *xc, const double *yc, const double *r)
        : InitialFunc() {SetParams(n, xc, yc, r);}

        virtual double XY(const double s, const double t) const;

    protected:

        virtual void SetParams(InputParams *params);
        void SetParams(const int n, const double *xc, const double *yc, const double *r);
    };

}

#endif
