/*************************************************
    um3boundary.h

    $Header: um3boundary.h,v 1.1 99/02/04 14:39:38 chopp Exp $

    $Log:       um3boundary.h,v $
    * Revision 1.1  99/02/04  14:39:38  14:39:38  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.1  97/12/04  10:30:02  10:30:02  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __UM3_BOUNDARY_H__
#define __UM3_BOUNDARY_H__

#include "boundary3d.h"
#include "uniformmesh3d.h"

namespace levelset {
        
    class UM3_Boundary : public Boundary3D
    {
    protected:
   
        UniformMesh3D* mesh;
   
    public:

    UM3_Boundary(const int wxlo = 1, const int wxhi = 1,
                 const int wylo = 1, const int wyhi = 1,
                 const int wzlo = 1, const int wzhi = 1)
        : Boundary3D(wxlo, wxhi, wylo, wyhi, wzlo, wzhi) {} 

        virtual void Apply(const int l) = 0;
        virtual void Apply(const int i, const int j, const int k, const int l) = 0;
        void SetMesh(UniformMesh3D* m) {mesh = m;}
        void SetIBoundary(const int l, const int n);
    };

}
#endif // __UM3_BOUNDARY_H__



