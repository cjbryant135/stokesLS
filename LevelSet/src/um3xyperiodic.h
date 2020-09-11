/*************************************************
    um3periodic.h

    $Header: um3periodic.h,v 1.1 99/02/04 14:39:43 chopp Exp $

    $Log:       um3periodic.h,v $
    * Revision 1.1  99/02/04  14:39:43  14:39:43  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  97/12/04  10:30:27  10:30:27  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __UM3XYPERIODIC_H__
#define __UM3XYPERIODIC_H__

#include "um3boundary.h"

namespace levelset {
        
    class UM3_XYPeriodicBdry : public UM3_Boundary 
    {
    public:

    UM3_XYPeriodicBdry(const int wxlo = 1, const int wxhi = 1,
                       const int wylo = 1, const int wyhi = 1,
                       const int wzlo = 1, const int wzhi = 1)
        : UM3_Boundary(wxlo, wxhi, wylo, wyhi, wzlo, wzhi) {} 
        virtual void Apply(const int l);
        virtual void Apply(const int i, const int j, const int k, const int l);
    };

}
#endif // __PERIODIC_H__
