/*************************************************
    um2boundary.h

    $Header: um2boundary.h,v 2.2 99/01/06 13:59:42 chopp Exp $

    $Log:       um2boundary.h,v $
    * Revision 2.2  99/01/06  13:59:42  13:59:42  chopp (David Chopp)
    * *** none ***
    * 
    * Revision 1.1  97/12/04  10:30:02  10:30:02  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __UM2_BOUNDARY_H__
#define __UM2_BOUNDARY_H__

#include "boundary2d.h"
#include "uniformmesh2d.h"

namespace levelset {
	 #ifndef __BCDIR_DEFINED__
	 #define __BCDIR_DEFINED__
    enum BCDir {Right=0, Left, Up, Down};
	 #endif

    class UM2_Boundary : public Boundary2D
    {
    protected:
   
        UniformMesh2D* mesh;
   
    public:

    UM2_Boundary(const int wxlo = 1, const int wxhi = 1,
                 const int wylo = 1, const int wyhi = 1)
        : Boundary2D(wxlo, wxhi, wylo, wyhi) {} 

        void SetMesh(UniformMesh2D* m) {mesh = m;}

        void SetIBoundary(const int k, const int n);
   
        int Nabor(const int i, const int j, const BCDir d, int& ni, int& nj);
    };

}

#endif // __BOUNDARY_H__
